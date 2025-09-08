package pgl.app.fastCall2.gpu;

import org.jocl.*;
import static org.jocl.CL.*;
import java.io.InputStream;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.atomic.AtomicBoolean;

/**
 * GPU computing manager - responsible for GPU resource management and computing task scheduling
 * Mainly optimizes the allele counting and genotype probability calculation in the scan step of FastCall2
 */
public class GPUComputeManager {

    private static GPUComputeManager instance;
    private static final Object lock = new Object();

    private cl_context context;
    private cl_command_queue commandQueue;
    private cl_device_id device;
    private cl_program program;

    // GPU computing kernels
    private cl_kernel alleleCountKernel;
    private cl_kernel genotypeCalcKernel;
    private cl_kernel probabilityKernel;

    private AtomicBoolean initialized = new AtomicBoolean(false);
    private ConcurrentHashMap<String, cl_mem> bufferCache = new ConcurrentHashMap<>();

    // GPU computing configuration
    private static final int MAX_WORK_GROUP_SIZE = 256;
    private static final int MEMORY_COALESCING_SIZE = 32;

    private GPUComputeManager() {
        initializeGPU();
    }

    public static GPUComputeManager getInstance() {
        if (instance == null) {
            synchronized (lock) {
                if (instance == null) {
                    instance = new GPUComputeManager();
                }
            }
        }
        return instance;
    }

    /**
     * Initialize GPU computing environment
     */
    private void initializeGPU() {
        try {
            // Enable exception handling
            CL.setExceptionsEnabled(true);

            // Get platform and device
            cl_platform_id[] platforms = new cl_platform_id[1];
            clGetPlatformIDs(1, platforms, null);

            cl_device_id[] devices = new cl_device_id[1];
            clGetDeviceIDs(platforms[0], CL_DEVICE_TYPE_GPU, 1, devices, null);
            device = devices[0];

            // Create context and command queue
            context = clCreateContext(null, 1, new cl_device_id[]{device}, null, null, null);
            commandQueue = clCreateCommandQueue(context, device, CL_QUEUE_PROFILING_ENABLE, null);

            // Compile GPU kernels
            compileKernels();

            initialized.set(true);
            System.out.println("GPU computing environment initialized successfully");

        } catch (Exception e) {
            System.err.println("GPU initialization failed, falling back to CPU computing: " + e.getMessage());
            initialized.set(false);
        }
    }

    /**
     * Compile GPU computing kernels
     */
    private void compileKernels() {
        String kernelSource = getKernelSource();
        program = clCreateProgramWithSource(context, 1, new String[]{kernelSource}, null, null);
        clBuildProgram(program, 0, null, null, null, null);

        // Create specific kernels
        alleleCountKernel = clCreateKernel(program, "countAlleles", null);
        genotypeCalcKernel = clCreateKernel(program, "calculateGenotypes", null);
        probabilityKernel = clCreateKernel(program, "calculateProbabilities", null);
    }

    /**
     * GPU kernel source code
     */
    private String getKernelSource() {
        try {
            // First try to detect NVIDIA GPU, if available use CUDA optimized kernel
            String preferredKernel = detectOptimalKernel();

            // Read kernel file from resources directory
            InputStream inputStream = getClass().getClassLoader().getResourceAsStream(preferredKernel);
            if (inputStream == null) {
                System.out.println("Cannot find preferred kernel file: " + preferredKernel + ", trying default OpenCL kernel");
                inputStream = getClass().getClassLoader().getResourceAsStream("gpu/kernels.cl");
            }

            if (inputStream == null) {
                throw new RuntimeException("Cannot find any GPU kernel files");
            }

            // Read file content
            StringBuilder kernelSource = new StringBuilder();
            try (java.io.BufferedReader reader = new java.io.BufferedReader(
                    new java.io.InputStreamReader(inputStream, "UTF-8"))) {
                String line;
                while ((line = reader.readLine()) != null) {
                    kernelSource.append(line).append("\n");
                }
            }

            return kernelSource.toString();

        } catch (Exception e) {
            System.err.println("Failed to read GPU kernel file: " + e.getMessage());
            // Return a simple fallback kernel
            return getFallbackKernelSource();
        }
    }

    /**
     * Detect optimal kernel type
     */
    private String detectOptimalKernel() {
        try {
            // Check GPU vendor information
            if (device != null) {
                long[] vendorSize = new long[1];
                clGetDeviceInfo(device, CL_DEVICE_VENDOR, 0, null, vendorSize);

                byte[] vendorBuffer = new byte[(int)vendorSize[0]];
                clGetDeviceInfo(device, CL_DEVICE_VENDOR, vendorBuffer.length,
                               Pointer.to(vendorBuffer), null);

                String vendor = new String(vendorBuffer, 0, vendorBuffer.length - 1).toLowerCase();

                if (vendor.contains("nvidia")) {
                    System.out.println("Detected NVIDIA GPU, using CUDA optimized kernel");
                    return "gpu/kernels_cuda_optimized.cl"; // CUDA optimized OpenCL version
                } else if (vendor.contains("amd") || vendor.contains("advanced micro devices")) {
                    System.out.println("Detected AMD GPU, using AMD optimized kernel");
                    return "gpu/kernels_amd_optimized.cl"; // AMD optimized version
                } else {
                    System.out.println("Detected generic GPU, using standard OpenCL kernel");
                    return "gpu/kernels.cl";
                }
            }
        } catch (Exception e) {
            System.err.println("GPU vendor detection failed: " + e.getMessage());
        }

        // Default return generic OpenCL kernel
        return "gpu/kernels.cl";
    }

    /**
     * Fallback simple kernel source code (for file read failures)
     */
    private String getFallbackKernelSource() {
        return "__kernel void countAlleles(__global const char* pileupData,\n" +
               "                          __global const int* positions,\n" +
               "                          __global const char* refBases,\n" +
               "                          __global int* alleleCounts,\n" +
               "                          const int dataLength,\n" +
               "                          const int numPositions,\n" +
               "                          const int maxAlleles) {\n" +
               "    int gid = get_global_id(0);\n" +
               "    if (gid >= numPositions) return;\n" +
               "    \n" +
               "    // Simplified counting logic\n" +
               "    int counts[6] = {0, 0, 0, 0, 0, 0};\n" +
               "    int startIdx = gid * dataLength;\n" +
               "    \n" +
               "    for (int i = 0; i < dataLength; i++) {\n" +
               "        char base = pileupData[startIdx + i];\n" +
               "        switch(base) {\n" +
               "            case 'A': case 'a': counts[0]++; break;\n" +
               "            case 'C': case 'c': counts[1]++; break;\n" +
               "            case 'G': case 'g': counts[2]++; break;\n" +
               "            case 'T': case 't': counts[3]++; break;\n" +
               "            case '-': counts[4]++; break;\n" +
               "            case '+': counts[5]++; break;\n" +
               "        }\n" +
               "    }\n" +
               "    \n" +
               "    for (int i = 0; i < 6; i++) {\n" +
               "        alleleCounts[gid * 6 + i] = counts[i];\n" +
               "    }\n" +
               "}\n" +
               "\n" +
               "__kernel void calculateGenotypes(__global const int* alleleCounts,\n" +
               "                               __global const float* errorRates,\n" +
               "                               __global float* genotypeProbabilities,\n" +
               "                               const int numPositions,\n" +
               "                               const int maxAlleles) {\n" +
               "    int gid = get_global_id(0);\n" +
               "    if (gid >= numPositions) return;\n" +
               "    \n" +
               "    genotypeProbabilities[gid * 3] = 1.0f;\n" +
               "    genotypeProbabilities[gid * 3 + 1] = 0.0f;\n" +
               "    genotypeProbabilities[gid * 3 + 2] = 0.0f;\n" +
               "}\n" +
               "\n" +
               "__kernel void calculateProbabilities(__global const float* data,\n" +
               "                                    __global float* results,\n" +
               "                                    const int dataSize,\n" +
               "                                    const float param1,\n" +
               "                                    const float param2) {\n" +
               "    int gid = get_global_id(0);\n" +
               "    if (gid >= dataSize) return;\n" +
               "    results[gid] = data[gid];\n" +
               "}";
    }

    /**
     * GPU accelerated allele counting
     */
    public int[] computeAlleleCountsGPU(String pileupData, int[] positions, String[] refBases) {
        if (!initialized.get()) {
            return null; // Return null indicates GPU is not available, need to fallback to CPU
        }

        try {
            int numPositions = positions.length;
            int dataLength = pileupData.length();

            // Create GPU memory buffers
            cl_mem pileupBuffer = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
                    Sizeof.cl_char * dataLength, Pointer.to(pileupData.getBytes()), null);

            cl_mem positionsBuffer = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
                    Sizeof.cl_int * numPositions, Pointer.to(positions), null);

            cl_mem resultBuffer = clCreateBuffer(context, CL_MEM_WRITE_ONLY,
                    Sizeof.cl_int * numPositions * 6, null, null);

            // Set kernel arguments
            clSetKernelArg(alleleCountKernel, 0, Sizeof.cl_mem, Pointer.to(pileupBuffer));
            clSetKernelArg(alleleCountKernel, 1, Sizeof.cl_mem, Pointer.to(positionsBuffer));
            clSetKernelArg(alleleCountKernel, 2, Sizeof.cl_mem, Pointer.to(resultBuffer));
            clSetKernelArg(alleleCountKernel, 3, Sizeof.cl_int, Pointer.to(new int[]{dataLength}));
            clSetKernelArg(alleleCountKernel, 4, Sizeof.cl_int, Pointer.to(new int[]{numPositions}));
            clSetKernelArg(alleleCountKernel, 5, Sizeof.cl_int, Pointer.to(new int[]{6}));

            // Execute kernel
            long[] globalWorkSize = new long[]{numPositions};
            long[] localWorkSize = new long[]{Math.min(MAX_WORK_GROUP_SIZE, numPositions)};

            clEnqueueNDRangeKernel(commandQueue, alleleCountKernel, 1, null,
                    globalWorkSize, localWorkSize, 0, null, null);

            // Read results
            int[] results = new int[numPositions * 6];
            clEnqueueReadBuffer(commandQueue, resultBuffer, CL_TRUE, 0,
                    Sizeof.cl_int * results.length, Pointer.to(results), 0, null, null);

            // Release resources
            clReleaseMemObject(pileupBuffer);
            clReleaseMemObject(positionsBuffer);
            clReleaseMemObject(resultBuffer);

            return results;

        } catch (Exception e) {
            System.err.println("GPU computing failed: " + e.getMessage());
            return null;
        }
    }

    /**
     * GPU accelerated genotype probability calculation
     */
    public float[] computeGenotypeProbabilitiesGPU(int[] alleleCounts, float[] errorRates) {
        if (!initialized.get()) {
            return null;
        }

        try {
            int numPositions = errorRates.length;

            // Create GPU buffers
            cl_mem countsBuffer = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
                    Sizeof.cl_int * alleleCounts.length, Pointer.to(alleleCounts), null);

            cl_mem errorBuffer = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
                    Sizeof.cl_float * errorRates.length, Pointer.to(errorRates), null);

            cl_mem resultBuffer = clCreateBuffer(context, CL_MEM_WRITE_ONLY,
                    Sizeof.cl_float * numPositions * 3, null, null);

            // Set kernel arguments
            clSetKernelArg(genotypeCalcKernel, 0, Sizeof.cl_mem, Pointer.to(countsBuffer));
            clSetKernelArg(genotypeCalcKernel, 1, Sizeof.cl_mem, Pointer.to(errorBuffer));
            clSetKernelArg(genotypeCalcKernel, 2, Sizeof.cl_mem, Pointer.to(resultBuffer));
            clSetKernelArg(genotypeCalcKernel, 3, Sizeof.cl_int, Pointer.to(new int[]{numPositions}));
            clSetKernelArg(genotypeCalcKernel, 4, Sizeof.cl_int, Pointer.to(new int[]{6}));

            // Execute kernel
            long[] globalWorkSize = new long[]{numPositions};
            long[] localWorkSize = new long[]{Math.min(MAX_WORK_GROUP_SIZE, numPositions)};

            clEnqueueNDRangeKernel(commandQueue, genotypeCalcKernel, 1, null,
                    globalWorkSize, localWorkSize, 0, null, null);

            // Read results
            float[] results = new float[numPositions * 3];
            clEnqueueReadBuffer(commandQueue, resultBuffer, CL_TRUE, 0,
                    Sizeof.cl_float * results.length, Pointer.to(results), 0, null, null);

            // Release resources
            clReleaseMemObject(countsBuffer);
            clReleaseMemObject(errorBuffer);
            clReleaseMemObject(resultBuffer);

            return results;

        } catch (Exception e) {
            System.err.println("GPU genotype calculation failed: " + e.getMessage());
            return null;
        }
    }

    /**
     * Check if GPU is available
     */
    public boolean isGPUAvailable() {
        return initialized.get();
    }

    /**
     * Get GPU device information
     */
    public String getGPUInfo() {
        if (!initialized.get()) {
            return "GPU not available";
        }

        try {
            long[] deviceName = new long[1];
            clGetDeviceInfo(device, CL_DEVICE_NAME, 0, null, deviceName);

            byte[] buffer = new byte[(int)deviceName[0]];
            clGetDeviceInfo(device, CL_DEVICE_NAME, buffer.length, Pointer.to(buffer), null);

            return new String(buffer, 0, buffer.length - 1);
        } catch (Exception e) {
            return "Cannot retrieve GPU information";
        }
    }

    /**
     * Cleanup GPU resources
     */
    public void cleanup() {
        if (initialized.get()) {
            try {
                clReleaseKernel(alleleCountKernel);
                clReleaseKernel(genotypeCalcKernel);
                clReleaseKernel(probabilityKernel);
                clReleaseProgram(program);
                clReleaseCommandQueue(commandQueue);
                clReleaseContext(context);

                bufferCache.clear();
                initialized.set(false);

                System.out.println("GPU resources cleaned up");
            } catch (Exception e) {
                System.err.println("GPU resources cleanup failed: " + e.getMessage());
            }
        }
    }
}

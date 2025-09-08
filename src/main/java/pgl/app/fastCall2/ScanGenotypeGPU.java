package pgl.app.fastCall2;

import pgl.app.fastCall2.gpu.GPUAlleleCounter;
import pgl.app.fastCall2.gpu.GPUComputeManager;
import pgl.infra.dna.FastaRecordBit;
import pgl.infra.utils.Benchmark;
import pgl.infra.utils.Dyad;
import pgl.infra.utils.PStringUtils;
import java.io.*;
import java.util.*;
import java.util.concurrent.*;
import java.util.concurrent.atomic.LongAdder;
import java.util.concurrent.atomic.AtomicLong;

/**
 * GPU optimization genotype scanning class
 * Mainly optimizes allele counting and genotype probability calculation in the scan step
 */
public class ScanGenotypeGPU extends ScanGenotype {

    private GPUAlleleCounter gpuCounter;
    private GPUComputeManager gpuManager;
    private boolean useGPUAcceleration;

    // Performance monitoring
    private AtomicLong gpuComputeTime = new AtomicLong(0);
    private AtomicLong cpuComputeTime = new AtomicLong(0);
    private AtomicLong totalPositionsProcessed = new AtomicLong(0);

    // GPU acceleration configuration
    private static final int GPU_BATCH_SIZE = 1000;
    private static final int MIN_DEPTH_FOR_GPU = 50;

    public ScanGenotypeGPU(String[] args) {
        super(args);
        initializeGPUComponents();
    }

    /**
     * Initialize GPU components
     */
    private void initializeGPUComponents() {
        try {
            this.gpuManager = GPUComputeManager.getInstance();
            this.gpuCounter = new GPUAlleleCounter();
            this.useGPUAcceleration = gpuManager.isGPUAvailable();

            if (useGPUAcceleration) {
                System.out.println("=== GPU heterogeneous optimization enabled ===");
                System.out.println("GPU Device: " + gpuManager.getGPUInfo());
                System.out.println("Batch size: " + GPU_BATCH_SIZE);
                System.out.println("GPU minimum depth threshold: " + MIN_DEPTH_FOR_GPU);
                System.out.println("========================");
            } else {
                System.out.println("GPU not available, using CPU computation mode");
            }
        } catch (Exception e) {
            System.err.println("GPU component initialization failed: " + e.getMessage());
            this.useGPUAcceleration = false;
        }
    }

    /**
     * GPU-optimized individual counting scan
     */
    @Override
    public void scanIndiCountsByThreadPool() {
        System.out.println("Starting GPU-optimized genotype scanning...");
        long startTime = System.nanoTime();

        // Using GPU-optimized IndiCount class
        FastaRecordBit frb = genomeFa.getFastaRecordBit(chromIndex);
        posRefMap = new HashMap<>();
        posAllelePackMap = new HashMap<>(vlEndIndex-vlStartIndex);
        positions = new int[vlEndIndex-vlStartIndex];

        for (int i = vlStartIndex; i < vlEndIndex; i++) {
            posRefMap.put(vl.positions[i], String.valueOf(frb.getBase(vl.positions[i]-1)));
            posAllelePackMap.put(vl.positions[i], vl.getAllelePacks(i));
            positions[i-vlStartIndex] = vl.positions[i];
        }

        Set<String> taxaSet = taxaBamsMap.keySet();
        ArrayList<String> taxaList = new ArrayList(taxaSet);
        Collections.sort(taxaList);

        Dyad<int[][], int[]> d = FastCall2.getBins(this.regionStart, this.regionEnd, FastCall2.scanBinSize);
        int[][] binBound = d.getFirstElement();
        int[] binStarts = d.getSecondElement();

        LongAdder counter = new LongAdder();
        ExecutorService pool = Executors.newFixedThreadPool(this.threadsNum);
        List<Future<IndiCount>> resultList = new ArrayList<>();

        for (int i = 0; i < taxaList.size(); i++) {
            List<String> bamPaths = taxaBamsMap.get(taxaList.get(i));
            StringBuilder sb = new StringBuilder(samtoolsPath);
            sb.append(" mpileup --no-output-ends ").append(this.baqMode).append("-q ").append(this.mappingQThresh).append(" -Q ").append(this.baseQThresh).append(" -f ").append(this.referenceFileS);
            for (int j = 0; j < bamPaths.size(); j++) {
                sb.append(" ").append(bamPaths.get(j));
            }
            sb.append(" -l ").append(vLibPosFileS).append(" -r ");
            sb.append(chrom).append(":").append(this.regionStart).append("-").append(this.regionEnd);
            String command = sb.toString();

            // Using GPU-optimized IndiCount
            IndiCountGPU idv = new IndiCountGPU(command, taxaList.get(i), binBound, binStarts, bamPaths, counter);
            Future<IndiCount> result = pool.submit(idv);
            resultList.add(result);
        }

        try {
            pool.shutdown();
            pool.awaitTermination(Long.MAX_VALUE, TimeUnit.SECONDS);
        }
        catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }

        long totalTime = System.nanoTime() - startTime;
        printPerformanceReport(totalTime);
    }

    /**
     * GPU optimization processing
     */
    private void optimizeWithGPU() {
        System.out.println("Executing GPU post-processing optimization...");

        try {
            // GPU-specific batch post-processing logic can be added here
            // For example: batch calculate genotype probabilities, optimize VCF generation, etc.

        } catch (Exception e) {
            System.err.println("GPU optimization processing failed: " + e.getMessage());
        }
    }

    /**
     * GPU-optimized IndiCount inner class
     */
    protected class IndiCountGPU extends IndiCount {

        private List<PendingComputation> pendingComputations;
        private int batchCounter = 0;

        public IndiCountGPU(String command, String taxonName, int[][] binBound,
                           int[] binStarts, List<String> bamPaths, LongAdder counter) {
            super(command, taxonName, binBound, binStarts, bamPaths, counter);
            this.pendingComputations = new ArrayList<>();
        }

        @Override
        public IndiCount call() throws Exception {
            System.out.println("Starting GPU-optimized individual counting - Sample: " + taxonName);
            long taskStartTime = System.nanoTime();

            try {
                Runtime rt = Runtime.getRuntime();
                Process p = rt.exec(command);

                BufferedReader br = new BufferedReader(new InputStreamReader(p.getInputStream()));
                BufferedReader bre = new BufferedReader(new InputStreamReader(p.getErrorStream()));

                String current = br.readLine();
                List<String> currentList = null;
                int currentPosition = -1;

                if (current != null) {
                    currentList = PStringUtils.fastSplit(current);
                    currentPosition = Integer.parseInt(currentList.get(1));
                }

                StringBuilder baseS = new StringBuilder();
                StringBuilder indelSb = new StringBuilder();

                // Process each position
                for (int i = 0; i < positions.length; i++) {
                    this.setDos(positions[i]);

                    if (current == null) {
                        this.writeMissing();
                    } else {
                        if (positions[i] == currentPosition) {
                            // Process current position
                            processPositionWithGPU(currentList, baseS, indelSb, currentPosition);

                            // Read next line
                            current = br.readLine();
                            if (current != null) {
                                currentList = PStringUtils.fastSplit(current);
                                currentPosition = Integer.parseInt(currentList.get(1));
                            }
                        } else if (positions[i] < currentPosition) {
                            this.writeMissing();
                        } else {
                            System.out.println("Position error, exiting program");
                            System.exit(1);
                        }
                    }
                }

                // Process remaining batch computations
                if (!pendingComputations.isEmpty()) {
                    processPendingComputations();
                }

                this.closeDos();
                br.close();
                p.waitFor();
                this.writeEmptyFiles();

                long taskTime = System.nanoTime() - taskStartTime;

                StringBuilder sb = new StringBuilder();
                sb.append("GPU-optimized individual counting completed - Sample: ").append(taxonName);
                sb.append(", Time taken: ").append(String.format("%.2f", taskTime / 1e9)).append(" seconds\n");

                // Add error messages
                String errorLine;
                while ((errorLine = bre.readLine()) != null) {
                    sb.append("Error: ").append(errorLine).append("\n");
                }

                System.out.println(sb.toString());

            } catch (Exception e) {
                System.err.println("Sample " + taxonName + " processing failed: " + e.getMessage());
                e.printStackTrace();
            }

            counter.increment();
            int cnt = counter.intValue();
            if (cnt % 10 == 0) {
                System.out.println(String.format("Completed %d/%d samples GPU-optimized counting",
                    cnt, taxaBamsMap.size()));
            }

            return this;
        }

        /**
         * Process position using GPU
         */
        private void processPositionWithGPU(List<String> currentList, StringBuilder baseS,
                                          StringBuilder indelSb, int currentPosition) {
            try {
                String ref = posRefMap.get(currentPosition);
                AllelePackage[] altAlleles = posAllelePackMap.get(currentPosition);

                baseS.setLength(0);
                int siteDepth = 0;

                // Merge data from all BAM files
                for (int j = 0; j < bamPaths.size(); j++) {
                    siteDepth += Integer.parseInt(currentList.get(3 + j * 3));
                    baseS.append(currentList.get(4 + j * 3));
                }

                if (siteDepth == 0) {
                    this.writeMissing();
                    return;
                }

                // Decide whether to use GPU based on depth
                if (useGPUAcceleration && siteDepth >= MIN_DEPTH_FOR_GPU) {
                    // Add to batch processing queue
                    pendingComputations.add(new PendingComputation(
                        altAlleles, baseS.toString(), siteDepth, currentPosition));

                    batchCounter++;

                    // Process when batch size is reached
                    if (batchCounter >= GPU_BATCH_SIZE) {
                        processPendingComputations();
                    }
                } else {
                    // Direct CPU processing
                    long cpuStart = System.nanoTime();
                    int[] alleleCounts = getAlleleCountsCPU(altAlleles, baseS.toString(), siteDepth, indelSb);
                    cpuComputeTime.addAndGet(System.nanoTime() - cpuStart);

                    this.writeAlleleCounts(alleleCounts);
                }

                totalPositionsProcessed.incrementAndGet();

            } catch (Exception e) {
                System.err.println("Position processing failed: " + e.getMessage());
                this.writeMissing();
            }
        }

        /**
         * Process pending GPU computations
         */
        private void processPendingComputations() {
            if (pendingComputations.isEmpty()) return;

            System.out.println("Executing GPU batch computation, batch size: " + pendingComputations.size());
            long gpuStart = System.nanoTime();

            try {
                // Prepare batch data
                List<GPUAlleleCounter.BatchInputData> batchData = new ArrayList<>();
                for (PendingComputation comp : pendingComputations) {
                    batchData.add(new GPUAlleleCounter.BatchInputData(
                        comp.altAlleles, comp.bases, comp.siteDepth));
                }

                // Execute GPU batch computation
                CompletableFuture<List<int[]>> future =
                    gpuCounter.batchProcessAlleleCountsAsync(batchData);

                List<int[]> results = future.get(30, TimeUnit.SECONDS); // 30 seconds timeout

                // Write results
                for (int i = 0; i < results.size() && i < pendingComputations.size(); i++) {
                    this.writeAlleleCounts(results.get(i));
                }

                gpuComputeTime.addAndGet(System.nanoTime() - gpuStart);

            } catch (Exception e) {
                System.err.println("GPU batch computation failed, switching to CPU: " + e.getMessage());

                // Fallback to CPU processing
                for (PendingComputation comp : pendingComputations) {
                    long cpuStart = System.nanoTime();
                    int[] alleleCounts = getAlleleCountsCPU(
                        comp.altAlleles, comp.bases, comp.siteDepth, new StringBuilder());
                    cpuComputeTime.addAndGet(System.nanoTime() - cpuStart);

                    this.writeAlleleCounts(alleleCounts);
                }
            } finally {
                pendingComputations.clear();
                batchCounter = 0;
            }
        }
    }

    /**
     * CPU version of allele counting (as fallback)
     */
    private int[] getAlleleCountsCPU(AllelePackage[] altAlleles, String bases,
                                   int siteDepth, StringBuilder indelSb) {
        // Call the CPU method of the GPU counter
        return gpuCounter.getAlleleCountsOptimized(altAlleles, bases, siteDepth, indelSb);
    }

    /**
     * Pending computation data structure
     */
    private static class PendingComputation {
        final AllelePackage[] altAlleles;
        final String bases;
        final int siteDepth;
        final int position;

        PendingComputation(AllelePackage[] altAlleles, String bases, int siteDepth, int position) {
            this.altAlleles = altAlleles;
            this.bases = bases;
            this.siteDepth = siteDepth;
            this.position = position;
        }
    }

    /**
     * Print performance report
     */
    private void printPerformanceReport(long totalTime) {
        System.out.println("\n=== GPU Optimization Performance Report ===");
        System.out.println("Total time: " + String.format("%.2f", totalTime / 1e9) + " seconds");
        System.out.println("GPU computation time: " + String.format("%.2f", gpuComputeTime.get() / 1e9) + " seconds");
        System.out.println("CPU computation time: " + String.format("%.2f", cpuComputeTime.get() / 1e9) + " seconds");
        System.out.println("Total positions processed: " + totalPositionsProcessed.get());

        if (gpuComputeTime.get() + cpuComputeTime.get() > 0) {
            double gpuRatio = (double)gpuComputeTime.get() / (gpuComputeTime.get() + cpuComputeTime.get());
            System.out.println("GPU computation ratio: " + String.format("%.1f%%", gpuRatio * 100));
        }

        if (useGPUAcceleration) {
            System.out.println("GPU acceleration status: " + gpuCounter.getPerformanceStats());
        }

        System.out.println("==================");
    }

    /**
     * Clean up resources
     */
    @Override
    protected void finalize() throws Throwable {
        if (gpuCounter != null) {
            gpuCounter.cleanup();
        }
        super.finalize();
    }
}

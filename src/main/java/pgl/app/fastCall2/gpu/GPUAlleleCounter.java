package pgl.app.fastCall2.gpu;

import pgl.app.fastCall2.AllelePackage;
import pgl.app.fastCall2.FastCall2;
import pgl.infra.utils.PStringUtils;
import java.util.concurrent.CompletableFuture;
import java.util.concurrent.ForkJoinPool;
import java.util.List;
import java.util.ArrayList;
import java.util.Arrays;

/**
 * GPU-accelerated allele counter
 * Optimizes the most time-consuming allele counting and genotype calculation steps in FastCall2
 */
public class GPUAlleleCounter {

    private final GPUComputeManager gpuManager;
    private final boolean useGPU;
    private final int batchSize;
    private static final int DEFAULT_BATCH_SIZE = 10000; // Batch processing size

    public GPUAlleleCounter() {
        this.gpuManager = GPUComputeManager.getInstance();
        this.useGPU = gpuManager.isGPUAvailable();
        this.batchSize = DEFAULT_BATCH_SIZE;

        if (useGPU) {
            System.out.println("GPU acceleration enabled - " + gpuManager.getGPUInfo());
        } else {
            System.out.println("GPU not available, using CPU computation mode");
        }
    }

    /**
     * GPU-optimized allele counting main method
     * This is the GPU-accelerated version of the original getAlleleCounts method
     */
    public int[] getAlleleCountsOptimized(AllelePackage[] altAlleles, String bases,
                                        int siteDepth, StringBuilder indelSb) {

        if (useGPU && siteDepth > 100) { // Use GPU only for sites with high depth
            return getAlleleCountsGPU(altAlleles, bases, siteDepth, indelSb);
        } else {
            return getAlleleCountsCPU(altAlleles, bases, siteDepth, indelSb);
        }
    }

    /**
     * GPU-accelerated allele counting
     */
    private int[] getAlleleCountsGPU(AllelePackage[] altAlleles, String bases,
                                   int siteDepth, StringBuilder indelSb) {
        try {
            // Preprocess data into GPU-friendly format
            GPUDataPreprocessor preprocessor = new GPUDataPreprocessor();
            GPUDataPreprocessor.GPUInputData inputData = preprocessor.prepareData(altAlleles, bases, siteDepth);

            // Use GPU to compute allele counts
            int[] gpuResults = gpuManager.computeAlleleCountsGPU(
                inputData.getPileupData(),
                inputData.getPositions(),
                inputData.getRefBases()
            );

            if (gpuResults != null) {
                // Post-process GPU results
                return preprocessor.processGPUResults(gpuResults, altAlleles, indelSb);
            } else {
                // GPU failed, fallback to CPU
                return getAlleleCountsCPU(altAlleles, bases, siteDepth, indelSb);
            }

        } catch (Exception e) {
            System.err.println("GPU computation error, switching to CPU: " + e.getMessage());
            return getAlleleCountsCPU(altAlleles, bases, siteDepth, indelSb);
        }
    }

    /**
     * CPU version of allele counting (original logic, used as fallback)
     */
    private int[] getAlleleCountsCPU(AllelePackage[] altAlleles, String bases,
                                   int siteDepth, StringBuilder indelSb) {

        // Preserve original CPU calculation logic to ensure normal operation when GPU is unavailable
        int[] alleleCounts = new int[altAlleles.length + 1]; // +1 for reference

        // Process simple substitutions
        for (int i = 0; i < bases.length(); i++) {
            char base = bases.charAt(i);

            // Reference allele (index 0)
            if (base == altAlleles[0].getRefAllele()) {
                alleleCounts[0]++;
                continue;
            }

            // Check alternative alleles
            boolean found = false;
            for (int j = 0; j < altAlleles.length; j++) {
                if (altAlleles[j].isSimpleAllele() &&
                    base == altAlleles[j].getAltAllele()) {
                    alleleCounts[j + 1]++;
                    found = true;
                    break;
                }
            }

            // Process Indels
            if (!found && (base == '+' || base == '-')) {
                processIndelCPU(bases, i, altAlleles, alleleCounts, indelSb);
            }
        }

        return alleleCounts;
    }

    /**
     * CPU version of Indel processing
     */
    private void processIndelCPU(String bases, int position, AllelePackage[] altAlleles,
                                int[] alleleCounts, StringBuilder indelSb) {

        indelSb.setLength(0);
        int i = position + 1;

        // Extract Indel sequence
        while (i < bases.length() && Character.isDigit(bases.charAt(i))) {
            indelSb.append(bases.charAt(i));
            i++;
        }

        if (indelSb.length() == 0) return;

        int indelLength = Integer.parseInt(indelSb.toString());
        indelSb.setLength(0);

        // Extract Indel sequence
        for (int j = 0; j < indelLength && i < bases.length(); j++, i++) {
            indelSb.append(bases.charAt(i));
        }

        String indelSeq = indelSb.toString().toUpperCase();

        // Match Indel alleles
        for (int j = 0; j < altAlleles.length; j++) {
            if (!altAlleles[j].isSimpleAllele()) {
                String altSeq = altAlleles[j].getAltAlleleString();
                if (indelSeq.equals(altSeq)) {
                    alleleCounts[j + 1]++;
                    break;
                }
            }
        }
    }

    /**
     * Batch processing of allele counts for multiple sites
     */
    public CompletableFuture<List<int[]>> batchProcessAlleleCountsAsync(
            List<BatchInputData> batchData) {

        return CompletableFuture.supplyAsync(() -> {
            List<int[]> results = new ArrayList<>();

            if (useGPU && batchData.size() > 10) {
                // GPU batch processing
                results = processBatchGPU(batchData);
            } else {
                // CPU parallel processing
                results = processBatchCPU(batchData);
            }

            return results;
        }, ForkJoinPool.commonPool());
    }

    /**
     * GPU batch processing
     */
    private List<int[]> processBatchGPU(List<BatchInputData> batchData) {
        List<int[]> results = new ArrayList<>();

        try {
            // Merge batch data into GPU-friendly format
            GPUBatchData gpuBatch = prepareBatchForGPU(batchData);

            // Execute GPU batch computation
            int[] batchResults = gpuManager.computeAlleleCountsGPU(
                gpuBatch.getCombinedPileupData(),
                gpuBatch.getPositionArray(),
                gpuBatch.getRefBaseArray()
            );

            if (batchResults != null) {
                // Parse batch results
                results = parseBatchGPUResults(batchResults, batchData);
            } else {
                // Fallback to CPU
                results = processBatchCPU(batchData);
            }

        } catch (Exception e) {
            System.err.println("GPU batch processing error, switching to CPU: " + e.getMessage());
            results = processBatchCPU(batchData);
        }

        return results;
    }

    /**
     * CPU parallel batch processing
     */
    private List<int[]> processBatchCPU(List<BatchInputData> batchData) {
        return batchData.parallelStream()
                .map(data -> getAlleleCountsCPU(
                    data.getAlleles(),
                    data.getBases(),
                    data.getDepth(),
                    new StringBuilder()))
                .toList();
    }

    /**
     * Prepare batch data for GPU
     */
    private GPUBatchData prepareBatchForGPU(List<BatchInputData> batchData) {
        StringBuilder combinedPileup = new StringBuilder();
        List<Integer> positions = new ArrayList<>();
        List<String> refBases = new ArrayList<>();

        for (int i = 0; i < batchData.size(); i++) {
            BatchInputData data = batchData.get(i);
            combinedPileup.append(data.getBases());
            positions.add(i); // Use index as position identifier

            // Get reference base
            if (data.getAlleles().length > 0) {
                refBases.add(String.valueOf(data.getAlleles()[0].getRefAllele()));
            } else {
                refBases.add("N");
            }
        }

        return new GPUBatchData(
            combinedPileup.toString(),
            positions.stream().mapToInt(Integer::intValue).toArray(),
            refBases.toArray(new String[0])
        );
    }

    /**
     * Parse GPU batch results
     */
    private List<int[]> parseBatchGPUResults(int[] batchResults, List<BatchInputData> batchData) {
        List<int[]> results = new ArrayList<>();
        int resultIndex = 0;

        for (BatchInputData data : batchData) {
            int numAlleles = data.getAlleles().length + 1;
            int[] siteResults = new int[numAlleles];

            // Extract results for this site from batch results
            for (int i = 0; i < Math.min(6, numAlleles); i++) {
                if (resultIndex < batchResults.length) {
                    siteResults[i] = batchResults[resultIndex++];
                }
            }

            results.add(siteResults);
        }

        return results;
    }

    /**
     * Batch input data structure
     */
    public static class BatchInputData {
        private final AllelePackage[] alleles;
        private final String bases;
        private final int depth;

        public BatchInputData(AllelePackage[] alleles, String bases, int depth) {
            this.alleles = alleles;
            this.bases = bases;
            this.depth = depth;
        }

        public AllelePackage[] getAlleles() { return alleles; }
        public String getBases() { return bases; }
        public int getDepth() { return depth; }
    }

    /**
     * GPU batch data structure
     */
    private static class GPUBatchData {
        private final String combinedPileupData;
        private final int[] positionArray;
        private final String[] refBaseArray;

        public GPUBatchData(String combinedPileupData, int[] positionArray, String[] refBaseArray) {
            this.combinedPileupData = combinedPileupData;
            this.positionArray = positionArray;
            this.refBaseArray = refBaseArray;
        }

        public String getCombinedPileupData() { return combinedPileupData; }
        public int[] getPositionArray() { return positionArray; }
        public String[] getRefBaseArray() { return refBaseArray; }
    }

    /**
     * Get performance statistics
     */
    public String getPerformanceStats() {
        if (useGPU) {
            return String.format("GPU acceleration mode - Device: %s", gpuManager.getGPUInfo());
        } else {
            return "CPU computation mode";
        }
    }

    /**
     * Clean up resources
     */
    public void cleanup() {
        if (useGPU) {
            gpuManager.cleanup();
        }
    }
}

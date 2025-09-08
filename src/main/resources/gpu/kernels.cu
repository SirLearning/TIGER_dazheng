// CUDA optimized kernels for FastCall2 GPU acceleration
// Optimized for NVIDIA GPUs using CUDA-specific features

#include <cuda_runtime.h>
#include <device_launch_parameters.h>

// CUDA optimized allele counting kernel
__global__ void countAllelesCUDA(const char* pileupData,
                                const int* positions,
                                const char* refBases,
                                int* alleleCounts,
                                const int dataLength,
                                const int numPositions,
                                const int maxAlleles) {

    int gid = blockIdx.x * blockDim.x + threadIdx.x;
    int lid = threadIdx.x;

    // Use shared memory for counting optimization (CUDA-specific)
    __shared__ int sharedCounts[256 * 6]; // 6 base types: A,C,G,T,-,+

    if (gid >= numPositions) return;

    // Initialize shared memory counters
    for (int i = lid; i < 256 * 6; i += blockDim.x) {
        sharedCounts[i] = 0;
    }
    __syncthreads();

    // Process pileup data
    int position = positions[gid];
    int startIdx = gid * dataLength;

    // Count various alleles
    for (int i = 0; i < dataLength; i++) {
        char base = pileupData[startIdx + i];
        int alleleType = -1;

        switch(base) {
            case 'A': case 'a': alleleType = 0; break;
            case 'C': case 'c': alleleType = 1; break;
            case 'G': case 'g': alleleType = 2; break;
            case 'T': case 't': alleleType = 3; break;
            case '-': alleleType = 4; break;
            case '+': alleleType = 5; break;
        }

        if (alleleType >= 0) {
            atomicAdd(&sharedCounts[lid * 6 + alleleType], 1);
        }
    }

    __syncthreads();

    // Merge shared counts to global memory
    if (lid == 0) {
        for (int i = 0; i < 6; i++) {
            int totalCount = 0;
            for (int j = 0; j < blockDim.x; j++) {
                totalCount += sharedCounts[j * 6 + i];
            }
            alleleCounts[gid * 6 + i] = totalCount;
        }
    }
}

// CUDA optimized genotype calculation kernel
__global__ void calculateGenotypesCUDA(const int* alleleCounts,
                                      const char* refBases,
                                      float* genotypeProbabilities,
                                      int* bestGenotypes,
                                      const int numPositions,
                                      const float qualityThreshold) {

    int gid = blockIdx.x * blockDim.x + threadIdx.x;
    if (gid >= numPositions) return;

    int baseIdx = gid * 6;

    // Get allele counts
    int countA = alleleCounts[baseIdx + 0];
    int countC = alleleCounts[baseIdx + 1];
    int countG = alleleCounts[baseIdx + 2];
    int countT = alleleCounts[baseIdx + 3];
    int countDel = alleleCounts[baseIdx + 4];
    int countIns = alleleCounts[baseIdx + 5];

    int totalDepth = countA + countC + countG + countT + countDel + countIns;

    if (totalDepth == 0) {
        bestGenotypes[gid] = -1; // Invalid genotype
        return;
    }

    // Calculate frequencies using CUDA math functions
    float freqA = __fdividef((float)countA, (float)totalDepth);
    float freqC = __fdividef((float)countC, (float)totalDepth);
    float freqG = __fdividef((float)countG, (float)totalDepth);
    float freqT = __fdividef((float)countT, (float)totalDepth);

    // Find the two highest frequency alleles
    float maxFreq1 = 0.0f, maxFreq2 = 0.0f;
    int allele1 = -1, allele2 = -1;

    float freqs[4] = {freqA, freqC, freqG, freqT};

    for (int i = 0; i < 4; i++) {
        if (freqs[i] > maxFreq1) {
            maxFreq2 = maxFreq1;
            allele2 = allele1;
            maxFreq1 = freqs[i];
            allele1 = i;
        } else if (freqs[i] > maxFreq2) {
            maxFreq2 = freqs[i];
            allele2 = i;
        }
    }

    // Calculate genotype probabilities (simplified likelihood calculation)
    float homozygousProb = maxFreq1 * maxFreq1;
    float heterozygousProb = 2.0f * maxFreq1 * maxFreq2;

    // Select the most likely genotype
    if (homozygousProb > heterozygousProb) {
        bestGenotypes[gid] = allele1 * 10 + allele1; // Homozygous (AA, CC, GG, TT)
        genotypeProbabilities[gid] = homozygousProb;
    } else {
        bestGenotypes[gid] = allele1 * 10 + allele2; // Heterozygous (AC, AG, AT, etc.)
        genotypeProbabilities[gid] = heterozygousProb;
    }
}

// CUDA optimized probability calculation kernel using Bayesian method
__global__ void calculateProbabilitiesCUDA(const int* alleleCounts,
                                          const float* qualityScores,
                                          float* posteriorProbabilities,
                                          const int numPositions,
                                          const float priorProb) {

    int gid = blockIdx.x * blockDim.x + threadIdx.x;
    if (gid >= numPositions) return;

    int baseIdx = gid * 6;

    // Get allele counts
    int counts[6];
    for (int i = 0; i < 6; i++) {
        counts[i] = alleleCounts[baseIdx + i];
    }

    int totalDepth = counts[0] + counts[1] + counts[2] + counts[3] + counts[4] + counts[5];

    if (totalDepth == 0) {
        posteriorProbabilities[gid] = 0.0f;
        return;
    }

    // Calculate likelihood function (binomial distribution approximation)
    float likelihood = 1.0f;
    float qualityScore = qualityScores ? qualityScores[gid] : 30.0f; // Default quality score
    float errorRate = __powf(10.0f, -qualityScore / 10.0f);

    // Likelihood calculation considering sequencing error rate
    for (int i = 0; i < 4; i++) { // Consider only ACGT
        if (counts[i] > 0) {
            float expectedErrorRate = __fdividef(errorRate, 3.0f); // Average distribution to other 3 bases
            float observedFreq = __fdividef((float)counts[i], (float)totalDepth);

            // Simplified likelihood calculation
            if (observedFreq > expectedErrorRate) {
                likelihood *= observedFreq;
            } else {
                likelihood *= expectedErrorRate;
            }
        }
    }

    // Bayesian posterior probability = likelihood × prior / evidence
    float posterior = likelihood * priorProb;

    // Normalization (simplified version)
    posteriorProbabilities[gid] = fminf(1.0f, posterior);
}

// CUDA optimized batch processing kernel
__global__ void batchProcessAllelesCUDA(const char* batchPileupData,
                                       const int* batchSizes,
                                       const int* batchOffsets,
                                       int* batchResults,
                                       const int numBatches) {

    int gid = blockIdx.x * blockDim.x + threadIdx.x;
    int batchIdx = blockIdx.y;

    if (batchIdx >= numBatches) return;

    int batchSize = batchSizes[batchIdx];
    int offset = batchOffsets[batchIdx];

    if (gid >= batchSize) return;

    // Process each position in the batch
    const char* pileupData = &batchPileupData[offset];
    int* results = &batchResults[batchIdx * batchSize * 6];

    // Execute allele counting
    int counts[6] = {0, 0, 0, 0, 0, 0};

    for (int i = 0; i < batchSize; i++) {
        char base = pileupData[gid * batchSize + i];
        switch(base) {
            case 'A': case 'a': counts[0]++; break;
            case 'C': case 'c': counts[1]++; break;
            case 'G': case 'g': counts[2]++; break;
            case 'T': case 't': counts[3]++; break;
            case '-': counts[4]++; break;
            case '+': counts[5]++; break;
        }
    }

    // Write results
    for (int i = 0; i < 6; i++) {
        results[gid * 6 + i] = counts[i];
    }
}

// CUDA optimized Indel detection and processing kernel
__global__ void processIndelsCUDA(const char* pileupData,
                                 const int* positions,
                                 int* indelInfo,
                                 char* indelSequences,
                                 const int dataLength,
                                 const int numPositions,
                                 const int maxIndelLength) {

    int gid = blockIdx.x * blockDim.x + threadIdx.x;
    if (gid >= numPositions) return;

    int startIdx = gid * dataLength;
    int indelCount = 0;
    int indelOutputIdx = gid * maxIndelLength;

    // Scan for indel markers in pileup data
    for (int i = 0; i < dataLength - 1; i++) {
        char current = pileupData[startIdx + i];

        if (current == '+' || current == '-') {
            // Parse indel length
            int indelLen = 0;
            int j = i + 1;

            // Read length digits
            while (j < dataLength && pileupData[startIdx + j] >= '0' && pileupData[startIdx + j] <= '9') {
                indelLen = indelLen * 10 + (pileupData[startIdx + j] - '0');
                j++;
            }

            // Record indel information
            if (indelLen > 0 && indelLen <= maxIndelLength && indelCount < 10) {
                indelInfo[gid * 10 + indelCount * 2] = (current == '+') ? 1 : -1; // Insertion: 1, Deletion: -1
                indelInfo[gid * 10 + indelCount * 2 + 1] = indelLen;

                // Copy indel sequence
                for (int k = 0; k < indelLen && j + k < dataLength; k++) {
                    indelSequences[indelOutputIdx + indelCount * maxIndelLength + k] = pileupData[startIdx + j + k];
                }

                indelCount++;
                i = j + indelLen - 1; // Skip processed sequence
            }
        }
    }

    // Fill remaining positions
    for (int i = indelCount; i < 10; i++) {
        indelInfo[gid * 10 + i * 2] = 0;
        indelInfo[gid * 10 + i * 2 + 1] = 0;
    }
}

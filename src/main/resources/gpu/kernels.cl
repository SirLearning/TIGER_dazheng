#pragma OPENCL EXTENSION cl_khr_global_int32_base_atomics : enable

// Allele counting kernel - core optimization algorithm
__kernel void countAlleles(__global const char* pileupData,
                         __global const int* positions,
                         __global const char* refBases,
                         __global int* alleleCounts,
                         const int dataLength,
                         const int numPositions,
                         const int maxAlleles) {

    int gid = get_global_id(0);
    int lid = get_local_id(0);
    int groupId = get_group_id(0);

    // Use local memory for counting optimization
    __local int localCounts[256 * 6]; // 6 base types: A,C,G,T,-,+

    if (gid >= numPositions) return;

    // Initialize local counters
    for (int i = lid; i < 256 * 6; i += get_local_size(0)) {
        localCounts[i] = 0;
    }
    barrier(CLK_LOCAL_MEM_FENCE);

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
            atomic_inc(&localCounts[lid * 6 + alleleType]);
        }
    }

    barrier(CLK_LOCAL_MEM_FENCE);

    // Merge local counts to global memory
    if (lid == 0) {
        for (int i = 0; i < 6; i++) {
            int totalCount = 0;
            for (int j = 0; j < get_local_size(0); j++) {
                totalCount += localCounts[j * 6 + i];
            }
            alleleCounts[gid * 6 + i] = totalCount;
        }
    }
}

// Genotype calculation kernel
__kernel void calculateGenotypes(__global const int* alleleCounts,
                               __global const char* refBases,
                               __global float* genotypeProbabilities,
                               __global int* bestGenotypes,
                               const int numPositions,
                               const float qualityThreshold) {

    int gid = get_global_id(0);
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

    // Calculate frequencies
    float freqA = (float)countA / totalDepth;
    float freqC = (float)countC / totalDepth;
    float freqG = (float)countG / totalDepth;
    float freqT = (float)countT / totalDepth;

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

// Probability calculation kernel - using Bayesian method
__kernel void calculateProbabilities(__global const int* alleleCounts,
                                   __global const float* qualityScores,
                                   __global float* posteriorProbabilities,
                                   const int numPositions,
                                   const float priorProb) {

    int gid = get_global_id(0);
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
    float errorRate = pow(10.0f, -qualityScore / 10.0f);

    // Likelihood calculation considering sequencing error rate
    for (int i = 0; i < 4; i++) { // Consider only ACGT
        if (counts[i] > 0) {
            float expectedErrorRate = errorRate / 3.0f; // Average distribution to other 3 bases
            float observedFreq = (float)counts[i] / totalDepth;

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
    posteriorProbabilities[gid] = fmin(1.0f, posterior);
}

// Indel detection and processing kernel
__kernel void processIndels(__global const char* pileupData,
                          __global const int* positions,
                          __global int* indelInfo,
                          __global char* indelSequences,
                          const int dataLength,
                          const int numPositions,
                          const int maxIndelLength) {

    int gid = get_global_id(0);
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

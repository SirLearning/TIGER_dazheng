package pgl.app.fastCall2.gpu;

import pgl.app.fastCall2.AllelePackage;
import pgl.app.fastCall2.FastCall2;
import java.util.HashMap;
import java.util.Map;
import java.util.regex.Pattern;
import java.util.regex.Matcher;

/**
 * GPU data preprocessor
 * Responsible for converting FastCall2's pileup data into GPU computation-friendly format
 */
public class GPUDataPreprocessor {

    // Base to numeric mapping (consistent with GPU kernel mapping)
    private static final Map<Character, Byte> BASE_TO_CODE = new HashMap<>();
    static {
        BASE_TO_CODE.put('A', (byte)0);
        BASE_TO_CODE.put('C', (byte)1);
        BASE_TO_CODE.put('G', (byte)2);
        BASE_TO_CODE.put('T', (byte)3);
        BASE_TO_CODE.put('-', (byte)4);
        BASE_TO_CODE.put('+', (byte)5);
        BASE_TO_CODE.put('a', (byte)0);
        BASE_TO_CODE.put('c', (byte)1);
        BASE_TO_CODE.put('g', (byte)2);
        BASE_TO_CODE.put('t', (byte)3);
    }

    // Indel pattern matching
    private static final Pattern INDEL_PATTERN = Pattern.compile("[+-](\\d+)([ACGTacgt]+)");

    /**
     * Prepare input data for GPU computation
     */
    public GPUInputData prepareData(AllelePackage[] altAlleles, String bases, int siteDepth) {
        // Preprocess pileup string, normalize format
        String normalizedBases = normalizePileupString(bases);

        // Create position array (for single site calculation, position is 0)
        int[] positions = {0};

        // Extract reference bases
        String[] refBases = new String[1];
        if (altAlleles.length > 0) {
            refBases[0] = String.valueOf(altAlleles[0].getRefAllele());
        } else {
            refBases[0] = "N";
        }

        return new GPUInputData(normalizedBases, positions, refBases);
    }

    /**
     * Normalize pileup string
     * Convert complex pileup format into GPU-friendly format for efficient processing
     */
    private String normalizePileupString(String bases) {
        if (bases == null || bases.isEmpty()) {
            return "";
        }

        StringBuilder normalized = new StringBuilder();
        int i = 0;

        while (i < bases.length()) {
            char c = bases.charAt(i);

            switch (c) {
                case '.':
                case ',':
                    // These symbols represent bases identical to the reference genome
                    // Will be specially marked during GPU processing
                    normalized.append('R'); // R represents Reference
                    i++;
                    break;

                case '^':
                    // Skip quality character
                    i += 2; // '^' + quality character
                    break;

                case '$':
                    // Read end marker, skip
                    i++;
                    break;

                case '*':
                    // Deletion marker
                    normalized.append('-');
                    i++;
                    break;

                case '+':
                case '-':
                    // Process Indel
                    String indelResult = processIndel(bases, i);
                    normalized.append(indelResult);
                    i = skipIndel(bases, i);
                    break;

                default:
                    // Regular bases A, C, G, T (case insensitive)
                    if (BASE_TO_CODE.containsKey(c)) {
                        normalized.append(Character.toUpperCase(c));
                    }
                    i++;
                    break;
            }
        }

        return normalized.toString();
    }

    /**
     * Process Indel sequences
     */
    private String processIndel(String bases, int startPos) {
        StringBuilder indelSeq = new StringBuilder();
        indelSeq.append(bases.charAt(startPos)); // '+' or '-'

        int i = startPos + 1;

        // Extract length digits
        while (i < bases.length() && Character.isDigit(bases.charAt(i))) {
            indelSeq.append(bases.charAt(i));
            i++;
        }

        // Extract sequence
        if (i < bases.length()) {
            try {
                int length = Integer.parseInt(indelSeq.substring(1));
                for (int j = 0; j < length && i < bases.length(); j++, i++) {
                    indelSeq.append(bases.charAt(i));
                }
            } catch (NumberFormatException e) {
                // Handle malformed indel, return marker only
                return String.valueOf(bases.charAt(startPos));
            }
        }

        return indelSeq.toString();
    }

    /**
     * Skip Indel sequence in original string
     */
    private int skipIndel(String bases, int startPos) {
        int i = startPos + 1;

        // Skip length digits
        while (i < bases.length() && Character.isDigit(bases.charAt(i))) {
            i++;
        }

        // Skip sequence
        if (i < bases.length()) {
            try {
                String lengthStr = bases.substring(startPos + 1, i);
                if (!lengthStr.isEmpty()) {
                    int length = Integer.parseInt(lengthStr);
                    i += length;
                }
            } catch (NumberFormatException e) {
                // Handle malformed indel
                i = startPos + 1;
            }
        }

        return i;
    }

    /**
     * Post-process GPU results
     */
    public int[] processGPUResults(int[] gpuResults, AllelePackage[] altAlleles, StringBuilder indelSb) {
        if (gpuResults == null || gpuResults.length == 0) {
            return new int[altAlleles.length + 1];
        }

        // GPU results format: [A_count, C_count, G_count, T_count, del_count, ins_count]
        int[] alleleCounts = new int[altAlleles.length + 1];

        // Map GPU base counts to allele counts
        for (int i = 0; i < altAlleles.length; i++) {
            AllelePackage allele = altAlleles[i];

            if (i == 0) {
                // Reference allele
                char refBase = allele.getRefAllele();
                alleleCounts[0] = getBaseCount(gpuResults, refBase);
            } else {
                // Alternative alleles
                if (allele.isSimpleAllele()) {
                    char altBase = allele.getAltAllele();
                    alleleCounts[i] = getBaseCount(gpuResults, altBase);
                } else {
                    // Complex alleles (Indels) handled separately
                    // For now, use deletion/insertion counts
                    String altSeq = allele.getAltAlleleString();
                    if (altSeq.startsWith("-")) {
                        alleleCounts[i] = gpuResults[4]; // Deletion count
                    } else if (altSeq.startsWith("+")) {
                        alleleCounts[i] = gpuResults[5]; // Insertion count
                    }
                }
            }
        }

        return alleleCounts;
    }

    /**
     * Get base count from GPU results
     */
    private int getBaseCount(int[] gpuResults, char base) {
        switch (Character.toUpperCase(base)) {
            case 'A': return gpuResults[0];
            case 'C': return gpuResults[1];
            case 'G': return gpuResults[2];
            case 'T': return gpuResults[3];
            default: return 0;
        }
    }

    /**
     * GPU input data container
     */
    public static class GPUInputData {
        private final String pileupData;
        private final int[] positions;
        private final String[] refBases;

        public GPUInputData(String pileupData, int[] positions, String[] refBases) {
            this.pileupData = pileupData;
            this.positions = positions;
            this.refBases = refBases;
        }

        public String getPileupData() { return pileupData; }
        public int[] getPositions() { return positions; }
        public String[] getRefBases() { return refBases; }
    }
}

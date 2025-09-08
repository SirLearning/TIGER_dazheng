package pgl.app.fastCall2;

import com.koloboke.collect.map.hash.HashByteByteMap;
import com.koloboke.collect.map.hash.HashByteByteMaps;
import pgl.AppUtils;
import pgl.infra.dna.allele.AlleleEncoder;
import pgl.infra.utils.*;
import java.util.*;

public class FastCall2 {

    String[] toolNames = {"disc","blib", "scan"};

    String currentTool = null;

    //genome block size for variation discovery. The max bin size should be less than 2^23. see {@link AllelePackage}
    static int disBinSize = 5000000;
    //genome block size for genotype scanning
    static int scanBinSize = 5000000;
    // maximum of the number of alternative alleles
    static int maxAltNum = 2;
    static int maxIndelLength = 63; // 就是一个long的长度
    //A, C, G, T, -, +
    static final byte[] pileupAlleleAscIIs = {65, 67, 71, 84, 45, 43};  // map to this 6, as {000, 001, 010, 011, 100, 101}

    static final HashByteByteMap pileupAscIIToAlleleCodingMap =
            HashByteByteMaps.getDefaultFactory().withDefaultValue((byte)-1).newImmutableMap(pileupAlleleAscIIs, AlleleEncoder.alleleCodings);   // what class withDefaultValue((byte)-1) returns?
    // key: pileupAlleleAscIIs; value: alleleCodings

    public FastCall2 (String[] args) {
        long timeStart = System.nanoTime();
        for (int i = 0; i < args.length; i++) {
            if (args[i].equals("-mod")) {
                this.currentTool = args[i+1];
                break;
            }
        }
        if (currentTool.equals("disc")) {
            System.out.println("Discovering genetic variation...");
            new DiscoverVariation(args);
        }
        else if (currentTool.equals("blib")) {
            System.out.println("Building the library of genetic variation from all samples...");
            new BuildVariationLibrary(args);
        }
        else if (currentTool.equals("vlib")) {
            System.out.println("View the library of genetic variation in text");
            new ViewVariationLibrary(args);
        }
        else if (currentTool.equals("clib")) {
            System.out.println("customize the library of genetic variation from position list");
            new CustomizeVariationLibrary(args);
        }
        else if (currentTool.equals("scan")) {
            System.out.println("Genotyping samples based on the variation library...");
            // Check if GPU optimization is enabled
            boolean useGPU = checkGPUOption(args);
            if (useGPU) {
                System.out.println("GPU heterogeneous optimization mode enabled");
                new ScanGenotypeGPU(args);
            } else {
                new ScanGenotype(args);
            }
        }
        else {
            System.out.println("Input errors in setting steps of FastCall 2. Programs stops.");
            System.exit(0);
        }
        StringBuilder sb = new StringBuilder("FastCall 2 is finished in ");
        sb.append((float)Benchmark.getTimeSpanHours(timeStart)).append(" hours.");
        System.out.println(sb.toString());
        System.out.println();
    }

    /**
     * Check if GPU optimization is enabled in command line arguments
     */
    private boolean checkGPUOption(String[] args) {
        for (int i = 0; i < args.length - 1; i++) {
            if (args[i].equals("-gpu") || args[i].equals("--gpu")) {
                String value = args[i + 1].toLowerCase();
                return value.equals("true") || value.equals("1") || value.equals("enable");
            }
        }
        // Enable GPU by default (if available)
        return true;
    }

    static Dyad<int[][], int[]> getBins (int regionStart, int regionEnd, int binSize) {
        int actualChrLength = regionEnd - regionStart;
        //starting from actual genome position
        int[][] binBound = PArrayUtils.getSubsetsIndicesBySubsetSize (actualChrLength, binSize);
        int[] binStarts = new int[binBound.length];
        for (int i = 0; i < binBound.length; i++) {
            binBound[i][0] = binBound[i][0]+regionStart;
            binBound[i][1] = binBound[i][1]+regionStart;
            binStarts[i] = binBound[i][0];
        }
        return new Dyad<>(binBound, binStarts);
    }

    static void removeFirstPositionSign (StringBuilder sb) {
        char charToRemove = '^';
        for (int i = 0; i < sb.length(); i++) {
            if (sb.charAt(i) == charToRemove) {
                sb.deleteCharAt(i);
                sb.deleteCharAt(i);
                i--;
            }
        }
    }
}

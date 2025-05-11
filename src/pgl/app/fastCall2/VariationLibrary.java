package pgl.app.fastCall2;

import it.unimi.dsi.fastutil.ints.IntArrayList;
import it.unimi.dsi.fastutil.ints.IntOpenHashSet;
import pgl.infra.utils.IOUtils;
import pgl.infra.utils.PArrayUtils;

import java.io.BufferedWriter;
import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.util.*;

public class VariationLibrary implements Comparable<VariationLibrary> {
    short chrom = Short.MIN_VALUE;
    //set to -1, when bin-based vls are concatenated to a chrom-based vl
    int binStart = Integer.MIN_VALUE;
    int[] positions = null;
    AllelePackage[][] allelePacks = null;
    private int[] potentialPositions = null;
    private ArrayList<AllelePackage>[] potentialAllelePackLists = null;


    public VariationLibrary (short chrom, int binStart, int[] positions, AllelePackage[][] allelePacks) {
        this.chrom = chrom;
        this.binStart = binStart;
        this.positions = positions;
        this.allelePacks = allelePacks;
    }

    public static VariationLibrary getInstance(List<VariationLibrary> vList) {
        Collections.sort(vList);    // VariationLibrary is comparable. compare binStart of instances.
        IntArrayList positionList = new IntArrayList();
        List<AllelePackage[]> alleleList = new ArrayList<>();
        for (int i = 0; i < vList.size(); i++) {
            for (int j = 0; j < vList.get(i).positions.length; j++) {
                positionList.add(vList.get(i).positions[j]);
                alleleList.add(vList.get(i).allelePacks[j]);
            }
        }
        int[] positions = positionList.toIntArray();    // combine positions in all VariationLibrary?
        // how to fix the bin problem? or it has been solved before?
        AllelePackage[][] allelePacks = alleleList.toArray(new AllelePackage[alleleList.size()][]); // combine allelePacks in all VariationLibrary
        VariationLibrary vl = new VariationLibrary(vList.get(0).chrom, -1, positions, allelePacks); // rearrange to a new VariationLibrary
        return vl;
    }

    public VariationLibrary (String infileS) {
        this.readBinaryFileS(infileS);
    }

    public VariationLibrary (List<IndividualGenotype> ingList, int maoThresh, int maxAltNum, short chrom, int binStart) {
        this.chrom = chrom;
        this.binStart = binStart;
        this.mergeIngs(ingList, maoThresh, maxAltNum);
    }

    private void mergeIngs (List<IndividualGenotype> ingList, int maoThresh, int maxAltNum) {
        IntOpenHashSet positionSet = new IntOpenHashSet();  // what is IntOpenHashSet? hashset is a set constructed in the hash way; but why it is open?
        for (int i = 0; i < ingList.size(); i++) {
            int positionNumber = ingList.get(i).getPositionNumber();    // how many allele positions this individual has
            for (int j = 0; j < positionNumber; j++) {
                positionSet.add(ingList.get(i).getAlleleChromPosition(j));  // the position is on the chromosome level, return an int
                // the set means no repeat value?
            }
        }
        potentialPositions = positionSet.toIntArray();  // allele positions, from a set (no repeat?) to int[]
        // the reason why these positions are potential is might that final allele positions will be filtered by threshold.
        Arrays.sort(potentialPositions);    // sort method for int[]
        potentialAllelePackLists = new ArrayList[potentialPositions.length];    // the arrayList's list number is based on position number
        for (int i = 0; i < potentialPositions.length; i++) {
            potentialAllelePackLists[i] = new ArrayList<>();    // why ArrayList? what is the use of ArrayList? is it the parent of every List class?
            // every potential position has a allelePackList
        }
        for (int i = 0; i < ingList.size(); i++) {
            int positionNum = ingList.get(i).getPositionNumber();   // how many allele positions this individual has
            int currentIndex = Integer.MIN_VALUE;
            for (int j = 0; j < positionNum; j++) {
                currentIndex = Arrays.binarySearch(potentialPositions, ingList.get(i).getAlleleChromPosition(j));   // get the index of matched positon in potentialPositions
                potentialAllelePackLists[currentIndex].add(new AllelePackage(ingList.get(i).getAllelePack(j)));
                // add the position-matched allelePack in the individual "i" to this allelePackList
            }
        }
        List<Integer> indexList = new ArrayList<>();    // the inner member of list must be class.
        AllelePackage[][] alts = new AllelePackage[potentialPositions.length][maxAltNum];
        // use AllelePackage as the data storage instance as it can store allelePack (int[]) and make list
        // dimension 1: every potential position has one;
        // dimension 2: maxAltNum = 2, means there are only 2 alternative allele at one position
        int[][] altCounts = new int[potentialPositions.length][maxAltNum];  // to count what? the individual (diploid) number of different alternative allele?
        // but is there a raw dataset to record maybe the all 4 allele? which is possible to happen
        for (int i = 0; i < potentialPositions.length; i++) {
            indexList.add(i);   // the reason why not use array like int[]? because array cannot expand length and cannot use methods in ArrayList.
        }
        indexList.parallelStream().forEach(i -> {   // learn more about the lambda method...
            Set<AllelePackage> alleleSet = new HashSet<>(); // the reason why AllelePackage implements Comparable? as set cannot have same members.
            for (int j = 0; j < potentialAllelePackLists[i].size(); j++) {
                alleleSet.add(potentialAllelePackLists[i].get(j));  // add allelePacks at this potential position (i) to alleleSet
            }
            AllelePackage[] alleles = alleleSet.toArray(new AllelePackage[alleleSet.size()]);   // why AllelePackage can also be used as array?
            Arrays.sort(alleles);   // array can use methods in Arrays class to do things like sort.
            int[] alleleCount = new int[alleles.length];    // why alleleCount have the length of AllelePackages number?
            int index = Integer.MIN_VALUE;
            for (int j = 0; j < potentialAllelePackLists[i].size(); j++) {
                index = Arrays.binarySearch(alleles, potentialAllelePackLists[i].get(j));   // find the same allelePack's index
                alleleCount[index]++;   // the int at the certain site of alleleCount represents the count number of same allele of alleles at the same site.
            }
            int[] countIndex = PArrayUtils.getIndicesByDescendingValue(alleleCount);    // get index order with alleleCount from big to small.
            int minNum = alleles.length;
            if (minNum > maxAltNum) minNum = maxAltNum; // there is only 2 alternative allele at maximum, also have more options (like 3, 4... in the future?)
            for (int j = 0; j < minNum; j++) {  // record to AllelePackage[][] and int[][]
                alts[i][j] = alleles[countIndex[j]];
                altCounts[i][j] = alleleCount[countIndex[j]];
            }
        });
        IntArrayList positionList = new IntArrayList(); // length is changeable.
        List<AllelePackage[]> allelePackLists = new ArrayList<AllelePackage[]>();   // to store allelePackList below?
        List<AllelePackage> allelePackList = new ArrayList<>();
        for (int i = 0; i < alts.length; i++) { // what the length is? position number?
            int varifiedNum = maxAltNum;    // means 2 at maximum.
            for (int j = maxAltNum; j > 0; j--) {
                if (altCounts[i][j-1] < maoThresh) varifiedNum--;   // check the minor allele occurrence (if < 2, abandon this allele)
            }
            if (varifiedNum == 0) continue;
            allelePackList.clear();
            for (int j = 0; j < varifiedNum; j++) {
                allelePackList.add(alts[i][j]); // add different alternative AllelePackage at the same potential site into a list.
            }
            positionList.add(alts[i][0].getAlleleChromPosition(binStart));  // every allelePackage record the position at "i".
            allelePackLists.add(allelePackList.toArray(new AllelePackage[allelePackList.size()]));
        }
        this.positions = positionList.toIntArray();
        this.allelePacks = allelePackLists.toArray(new AllelePackage[allelePackLists.size()][]);
        // why not record numbers for each alternative allele?
    }

    public AllelePackage[] getAllelePacks (int positionIndex) {
        return this.allelePacks[positionIndex];
    }

    public int getPositionIndex (int pos) {
        return Arrays.binarySearch(positions, pos);
    }

    /**
     * Return the index of the next position in the library, inclusive
     * @param pos
     * @return Integer.MIN_VALUE if the query pos is greater than the end of the position list
     */
    public int getStartIndex (int pos) {
        int index = this.getPositionIndex(pos);
        if (index < 0) {
            index = -index -1;
            if (index == positions.length) {
                return Integer.MIN_VALUE;
            }
        }
        return index;
    }

    /**
     * Return the index of the previous position in the library, exclusive
     * @param pos
     * @return Integer.MIN_VALUE if the query pos is less than the start of the position list
     */
    public int getEndIndex (int pos) {
        int index = this.getPositionIndex(pos);
        if (index < 0) {
            index = -index -2;
            if (index < 0) return Integer.MIN_VALUE;
        }
        return index+1;
    }

    public short getChrom () {
        return this.chrom;
    }

    public int getPosition (int positionIndex) {
        return this.positions[positionIndex];
    }

    public void writeBinaryFileS (String outfileS) {
        try {
            DataOutputStream dos = IOUtils.getBinaryGzipWriter(outfileS);
            dos.writeShort(chrom);
            dos.writeInt(binStart); // why start with -1?
            dos.writeInt(positions.length);
            for (int i = 0; i < positions.length; i++) {
                dos.writeInt(positions[i]);
                dos.writeByte((byte) this.allelePacks[i].length);   //  < maxAltNum
                for (int j = 0; j < allelePacks[i].length; j++) {
                    int num = allelePacks[i][j].getAllelePackSize();    // will be different when the allelePack is an indel
                    for (int k = 0; k < num; k++) {
                        dos.writeInt(allelePacks[i][j].getAllelePack()[k]); // write int one by one
                        // why not record occurrence number of every allelePack? is it will be counted after again?
                    }
                }
            }
            dos.flush();
            dos.close();
        }
        catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
        StringBuilder sb = new StringBuilder();
        sb.append(positions.length).append(" polymorphic sites are written to ").append(outfileS);
        System.out.println(sb.toString());
    }

    public void writeBinaryFileS (String outfileS, int[] customPositions) {
        IntArrayList indexList = new IntArrayList();
        for (int i = 0; i < positions.length; i++) {
            int index = Arrays.binarySearch(customPositions, positions[i]);
            if (index < 0) continue;
            indexList.add(i);
        }
        int[] indices = indexList.toIntArray();
        try {
            DataOutputStream dos = IOUtils.getBinaryGzipWriter(outfileS);
            dos.writeShort(chrom);
            dos.writeInt(binStart);
            dos.writeInt(indices.length);
            for (int i = 0; i < indices.length; i++) {
                dos.writeInt(positions[indices[i]]);
                dos.writeByte((byte) this.allelePacks[indices[i]].length);
                for (int j = 0; j < allelePacks[indices[i]].length; j++) {
                    int num = allelePacks[indices[i]][j].getAllelePackSize();
                    for (int k = 0; k < num; k++) {
                        dos.writeInt(allelePacks[indices[i]][j].getAllelePack()[k]);
                    }
                }
            }
            dos.flush();
            dos.close();
        }
        catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
        StringBuilder sb = new StringBuilder();
        sb.append(positions.length).append(" polymorphic sites are written to ").append(outfileS);
        System.out.println(sb.toString());
    }

    private void readBinaryFileS (String infileS) {
        try {
            DataInputStream dis = IOUtils.getBinaryGzipReader(infileS);
            chrom = dis.readShort();
            binStart = dis.readInt();
            int positionNum = dis.readInt();
            positions = new int[positionNum];
            allelePacks = new AllelePackage[positionNum][];
            for (int i = 0; i < positionNum; i++) {
                positions[i] = dis.readInt();
                int alleleNum = dis.readByte();
                allelePacks[i] = new AllelePackage[alleleNum];
                for (int j = 0; j < alleleNum; j++) {
                    int firstInt = dis.readInt();
                    int packSize = AllelePackage.getAllelePackSizeFromFirstInt(firstInt);
                    int[] allelePack = new int[packSize];
                    allelePack[0] = firstInt;
                    for (int k = 1; k < packSize; k++) {
                        allelePack[k] = dis.readInt();
                    }
                    allelePacks[i][j] = new AllelePackage(allelePack);
                }
            }
        }
        catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
    }

    public void writeTextFileS (String outfileS) {
        try {
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            bw.write("Position\tAlts\tIndelLength\tIndelSeq");
            bw.newLine();
            StringBuilder sb = new StringBuilder();
            for (int i = 0; i < positions.length; i++) {
                sb.setLength(0);
                sb.append(positions[i]).append("\t");
                for (int j = 0; j < this.allelePacks[i].length; j++) {
                    sb.append(allelePacks[i][j].getAlleleBase()).append(",");
                }
                sb.deleteCharAt(sb.length()-1).append("\t");
                int totalLength = 0;
                for (int j = 0; j < this.allelePacks[i].length; j++) {
                    int indelLength = allelePacks[i][j].getIndelLength();
                    sb.append(indelLength).append(",");
                    totalLength+=indelLength;
                }
                sb.deleteCharAt(sb.length()-1);
                if (totalLength > 0) {
                    sb.append("\t");
                    for (int j = 0; j < this.allelePacks[i].length; j++) {
                        sb.append(allelePacks[i][j].getIndelSeq()).append(",");
                    }
                    sb.deleteCharAt(sb.length()-1);
                }
                bw.write(sb.toString());
                bw.newLine();
            }
            bw.flush();
            bw.newLine();
        }
        catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
    }

    @Override
    public int compareTo(VariationLibrary o) {
        if (this.binStart < o.binStart) return -1;
        else if (this.binStart > o.binStart) return 1;
        return 0;
    }
}

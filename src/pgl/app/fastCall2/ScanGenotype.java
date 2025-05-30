package pgl.app.fastCall2;

import com.koloboke.collect.map.IntDoubleMap;
import com.koloboke.collect.map.hash.HashIntDoubleMaps;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.io.FileUtils;
import pgl.AppAbstract;
import pgl.PGLAPPEntrance;
import pgl.PGLConstraints;
import pgl.infra.dna.FastaBit;
import pgl.infra.dna.FastaRecordBit;
import pgl.infra.dna.allele.AlleleEncoder;
import pgl.infra.utils.Benchmark;
import pgl.infra.utils.Dyad;
import pgl.infra.utils.IOUtils;
import pgl.infra.utils.PStringUtils;
import java.io.*;
import java.text.SimpleDateFormat;
import java.util.*;
import java.util.concurrent.*;
import java.util.concurrent.atomic.LongAdder;

import static cern.jet.math.Arithmetic.factorial;

class ScanGenotype extends AppAbstract {
    //Reference genome file with an index file (.fai). The reference should be in Fasta format. Chromosomes are labled as 1-based numbers (1,2,3,4,5...).
    String referenceFileS = null;
    //The taxaBamMap file contains information of taxon and its corresponding bam files. The bam file should have .bai file in the same folder.
    String taxaBamMapFileS = null;
    //The genetic variation library file
    String libFileS = null;
    int chrom = Integer.MIN_VALUE;
    //Starting position of the specified region for variation calling, inclusive
    int regionStart = Integer.MIN_VALUE;
    //Ending position the specified regionfor variation calling, exclusive
    int regionEnd = Integer.MIN_VALUE;
    //The switch of base alignment quality (BAQ) computaiton, 0 is diabled and 1 is enbabled.
    String baqMode = "-B ";
    //Minimum mapping quality (MQ) for an alignment to be used for genotyping.
    int mappingQThresh = 30;
    //Minimum base quality (BQ) for a base to be used for genotyping.
    int baseQThresh = 20;
    //combined: sequencing error and alignment error
    double combinedErrorRate = 0.05;
    //The path of samtools
    String samtoolsPath = null;
    //VCF output directory
    String outputDirS = null;
    //Number of threads (taxa number to be processed at the same time)
    int threadsNum = PGLConstraints.parallelLevel;

    String[] subDirS = {"indiVCF", "indiCounts", "VCF"};    // why there is indiVCF?

    HashMap<String, List<String>> taxaBamsMap = null;
    HashMap<String, Double> taxaCoverageMap = null;
    String[] taxaNames = null;

    IntDoubleMap factorialMap = null;
    int maxFactorial = 150; // why 150? and is used for what? short reads alignment? 150 is not a factorial?

    String vLibPosFileS = null;

    FastaBit genomeFa = null;
    int chromIndex = Integer.MIN_VALUE;
    VariationLibrary vl = null;
    int vlStartIndex = Integer.MIN_VALUE;
    int vlEndIndex = Integer.MIN_VALUE;
    HashMap<Integer, String> posRefMap = new HashMap<>();
    //HashMap would be faster to locate alleles, larger memory though
    HashMap<Integer, AllelePackage[]> posAllelePackMap = new HashMap<>();
    int[] positions = null;
    int vlBinStartIndex = 0;
    int vlBinEndIndex = 0;

    public ScanGenotype (String[] args) {
        this.creatAppOptions();
        this.retrieveAppParameters(args);

        this.mkDir();   // make 3 directions
        this.processVariationLibrary(); // create a file store positions in the region
        this.creatFactorialMap();

        /*
        Output by individual allele count, fast
         */
        this.scanIndiCountsByThreadPool();
        this.mkFinalVCFFromIndiCounts();
    }

    @Override
    public void creatAppOptions() {
        options.addOption("app", true, "App name.");
        options.addOption("mod", true, "Module name of FastCall 2.");
        options.addOption("a", true, "Reference genome file with an index file (.fai). The reference should be in Fasta format. " +
            "Chromosomes are labelled as numbers (1,2,3,4,5...). It is recommended to use reference chromosome while perform genotyping for " +
            "each chromosome because loading reference genome would be much faster.");
        options.addOption("b", true, "The taxaBamMap file contains information of taxon and its corresponding bam files. " +
            "The bam file should have .bai file in the same folder.");
        options.addOption("c", true, "The genetic variation library file.");
        options.addOption("d", true, "Chromosome or region on which genotyping will be performed (e.g. chromosome 1 is designated as 1. " +
            "Region 1bp to 100000bp on chromosome 1 is 1:1,100000)");
        options.addOption("e", true, "The switch of base alignment quality (BAQ) computaiton, 0 is diabled and 1 is enbabled. It is 0 by default.");
        options.addOption("f", true, "Minimum mapping quality (MQ) for an alignment to be used for genotyping. It is 30 by default.");
        options.addOption("g", true, "Minimum base quality (BQ) for a base to be used for genotyping. It is 20 by default.");
        options.addOption("h", true, "Combined error rate of sequencing and misalignment. Heterozygous read mapping are more " +
            "likely to be genotyped as homozygote when the combined error rate is high.");
        options.addOption("i", true, "The path of samtools.");
        options.addOption("j", true, "Number of threads. It is 32 by default.");
        options.addOption("k", true, "The directory of VCF output.");
    }

    @Override
    public void retrieveAppParameters(String[] args) {
        CommandLineParser parser = new DefaultParser();
        try {
            CommandLine line = parser.parse(options, args);
            String inOpt = null;
            this.referenceFileS = line.getOptionValue("a");
            this.taxaBamMapFileS = line.getOptionValue("b");
            this.libFileS = line.getOptionValue("c");
            String[] tem = line.getOptionValue("d").split(":");
            this.chrom = Integer.parseInt(tem[0]);
            long start = System.nanoTime();
            System.out.println("Reading reference genome from "+ referenceFileS);
            this.genomeFa = new FastaBit(referenceFileS);
            System.out.println("Reading reference genome took " + String.format("%.2f", Benchmark.getTimeSpanSeconds(start)) + "s");
            this.chromIndex = genomeFa.getIndexByDescription(String.valueOf(this.chrom));
            if (tem.length == 1) {
                this.regionStart = 1;
                this.regionEnd = genomeFa.getSeqLength(chromIndex)+1;
            }
            else if (tem.length == 2) {
                tem = tem[1].split(",");
                this.regionStart = Integer.parseInt(tem[0]);
                this.regionEnd = Integer.parseInt(tem[1])+1;
            }
            inOpt = line.getOptionValue("e");
            if (inOpt != null) {
                if (inOpt.equals("1")) {
                    this.baqMode ="";
                }
                else if (inOpt.equals("0")) {
                    this.baqMode = "-B ";
                }
                else {
                    System.out.println("Warning: Incorrect input for -c option. The BAQ computation is disabled by default");
                }
                inOpt = null;
            }
            inOpt = line.getOptionValue("f");
            if (inOpt != null) {
                this.mappingQThresh = Integer.parseInt(inOpt);
                inOpt = null;
            }
            inOpt = line.getOptionValue("g");
            if (inOpt != null) {
                this.baseQThresh = Integer.parseInt(inOpt);
                inOpt = null;
            }
            inOpt = line.getOptionValue("h");
            if (inOpt != null) {
                this.combinedErrorRate = Double.parseDouble(inOpt);
                inOpt = null;
            }
            this.samtoolsPath = line.getOptionValue("i");
            inOpt = line.getOptionValue("j");
            if (inOpt != null) {
                this.threadsNum = Integer.parseInt(inOpt);
                inOpt = null;
            }
            this.outputDirS = line.getOptionValue("k");
        }
        catch(Exception e) {
            e.printStackTrace();
            System.out.println("\nThere are input errors in the command line. Program stops.");
            this.printInstructionAndUsage();
            System.exit(0);
        }
        this.parseTaxaBamMap(this.taxaBamMapFileS); // parse the taxa and bam information
        // the taxon name, coverage and bam file sites are recorded in the taxaBamMap.txt file
    }

    @Override
    public void printInstructionAndUsage() {
        System.out.println(PGLAPPEntrance.getTIGERIntroduction());
        System.out.println("Below are the commands of module \"scan\" in FastCall 2.");
        this.printUsage();
    }

    public void mkFinalVCFFromIndiCounts () {
        String outfileS = new File(outputDirS, subDirS[2]).getAbsolutePath();
        outfileS = new File(outfileS, "chr"+PStringUtils.getNDigitNumber(3, chrom)+".vcf").getAbsolutePath();   // finally!
        try {
            SimpleDateFormat sdf = new SimpleDateFormat("MM/dd/yyyy HH:mm:ss.SSS"); // format method to record the date
            Date dt = new Date();
            String S = sdf.format(dt);
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            bw.write("##fileformat=VCFv4.1\n");
            bw.write("##fileDate="+S.split(" ")[0]+"\n");
            bw.write("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n");
            bw.write("##FORMAT=<ID=AD,Number=.,Type=Integer,Description=\"Allelic depths for the reference and alternate alleles in the order listed\">\n");
            bw.write("##FORMAT=<ID=GL,Number=G,Type=Integer,Description=\"Genotype likelihoods for 0/0, 0/1, 1/1, or  0/0, 0/1, 0/2, 1/1, 1/2, 2/2 if 2 alt alleles\">\n");
            bw.write("##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">\n");
            bw.write("##INFO=<ID=NZ,Number=1,Type=Integer,Description=\"Number of taxa with called genotypes\">\n");
            bw.write("##INFO=<ID=AD,Number=.,Type=Integer,Description=\"Total allelelic depths in order listed starting with REF\">\n");
            bw.write("##INFO=<ID=AC,Number=.,Type=Integer,Description=\"Numbers of allele occurence across taxa in order listed\">\n");
            bw.write("##INFO=<ID=IS,Number=.,Type=Integer,Description=\"Indel sequence of ALT alleles in order listed\">\n");
            bw.write("##INFO=<ID=GN,Number=.,Type=Integer,Description=\"Number of taxa with genotypes AA,AB,BB or AA,AB,AC,BB,BC,CC if 2 alt alleles\">\n");
            bw.write("##INFO=<ID=HT,Number=1,Type=Integer,Description=\"Number of heterozygotes\">\n");
            bw.write("##INFO=<ID=MAF,Number=1,Type=Float,Description=\"Minor allele frequency\">\n");
            bw.write("##ALT=<ID=D,Description=\"Deletion\">\n");
            bw.write("##ALT=<ID=I,Description=\"Insertion\">\n");
            Dyad<int[][], int[]> d = FastCall2.getBins(this.regionStart, this.regionEnd, FastCall2.scanBinSize);    // scanBinSize is also 5M bp
            int[][] binBound = d.getFirstElement();
            int[] binStarts = d.getSecondElement(); // used to identify the number of bin
                StringBuilder sb = new StringBuilder("#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT");
            for (int i = 0; i < taxaNames.length; i++) {
                sb.append("\t").append(taxaNames[i]);
            }
            bw.write(sb.toString());
            bw.newLine();
            List<Future<IndividualCount>> futureList = new ArrayList<>();
            List<IndividualCount> incList = new ArrayList<>();  // store IndividualCount for each taxon at each bin
            vlBinStartIndex = 0;
            vlBinEndIndex = 0;
            for (int i = 0; i < binStarts.length; i++) {    // recursive for every bin
                futureList.clear();
                incList.clear();
                String indiCountFolderS = new File(outputDirS, subDirS[1]).getAbsolutePath();
                try {
                    LongAdder counter = new LongAdder();
                    ExecutorService pool = Executors.newFixedThreadPool(this.threadsNum);
                    sb.setLength(0);
                    sb.append(chrom).append("_").append(binBound[i][0]).append("_").append(binBound[i][1]).append(".iac.gz");   // read .iac in this bin
                    // why read again? not directly store in memory?
                    for (int j = 0; j < taxaNames.length; j++) {    // do for every taxon
                        String indiTaxonDirS = new File (indiCountFolderS, taxaNames[j]).getAbsolutePath();
                        String fileS = new File (indiTaxonDirS, sb.toString()).getAbsolutePath();
                        TaxonCountRead tr = new TaxonCountRead(fileS);  // read a .iac file for each taxon and each bin
                        Future<IndividualCount> result = pool.submit(tr);   // only to store IndividualCount into memory, transfer from one method to this one
                        // IndividualCount is not same as IndiCount
                        futureList.add(result);
                    }
                    pool.shutdown();
                    pool.awaitTermination(Long.MAX_VALUE, TimeUnit.SECONDS);    // wait util all read finish
                    for (int j = 0; j < futureList.size(); j++) {
                        IndividualCount inc = futureList.get(j).get();
                        if (inc == null) continue;  // IndividualCount implements Comparable
                        incList.add(inc);
                    }
                    Collections.sort(incList);  // sort the taxa name
                }
                catch (Exception e) {
                    e.printStackTrace();
                    System.exit(1);
                }
                vlBinEndIndex = vlBinStartIndex + incList.get(0).alleleNum.length;  // position number as index number of VariationLibrary
                List<Integer> indexList = new ArrayList<>();
                for (int j = 0; j < incList.get(0).alleleNum.length; j++) {
                    indexList.add(j);   // a list of the relative positon index (from 0) in the bin
                }
                String[] vcfRecords = new String[indexList.size()]; // record the whole VCF file
                indexList.parallelStream().forEach(index -> {   // parallel for each position
                    StringBuilder vsb = new StringBuilder();    // lines in VCF file
                    int currentPosition = positions[index+vlBinStartIndex]; // from relative index in bin to real position on chromosome
                    // write information of this position (in 1 line)
                    vsb.append(chrom).append("\t").append(currentPosition).append("\t").append(chrom).append("-").append(currentPosition)
                            .append("\t").append(posRefMap.get(currentPosition)).append("\t");
                    // begin with: CHROM POS ID REF
                    AllelePackage[] altAlleles = posAllelePackMap.get(currentPosition);
                    for (int j = 0; j < altAlleles.length; j++) {
                        vsb.append(altAlleles[j].getAlleleBase()).append(",");  // add ALT (seperated by ",")
                    }
                    vsb.deleteCharAt(vsb.length()-1).append("\t.\t.\t");    // delete "," and add QUAL and FILTER
                    List<short[]> siteCountsList = new ArrayList<>();
                    for (int j = 0; j < incList.size(); j++) {  // record for each taxon
                        siteCountsList.add(incList.get(j).alleleCounts[index]); // get alleleCounts[index], is a length changeable short[] (length is according to alternative allele number, including reference genotype)
                    }
                    vsb.append(this.getInfoAndGenotypes(siteCountsList, altAlleles));   // format the INFO and Genotype elements in the line
                    vcfRecords[index] = vsb.toString();
                });
                for (int j = 0; j < vcfRecords.length; j++) {
                    bw.write(vcfRecords[j]);    // write lines
                    bw.newLine();
                }
                vlBinStartIndex = vlBinEndIndex;
                System.out.println(sb.toString().split("\\.")[0]+ " genotyping is finished.");
            }
            bw.flush();
            bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
        this.deleteTemperateFile(); // what files?
        System.out.println("Final VCF is completed at " + outfileS);
        System.out.println("Genotyping is finished.");
    }

    class TaxonCountRead implements Callable<IndividualCount> {
        String fileS;
        public TaxonCountRead (String fileS) {
            this.fileS = fileS;
        }

        @Override
        public IndividualCount call() throws Exception {
            File f = new File (fileS);
            if (!f.exists()) {
                System.out.println("Warning: "+ f.getAbsolutePath()+" does not exist");
                return null;
            }
            IndividualCount inc = new IndividualCount(this.fileS);  // readBinaryFileS
            // maybe the time to read file from hardware is long. therefore, use parallel
            return inc;
        }
    }

    public void scanIndiCountsByThreadPool () {
        FastaRecordBit frb = genomeFa.getFastaRecordBit(chromIndex);    // store by Sequence3Bit field, which is inherited by FastaRecordBit
        posRefMap = new HashMap<>();    // reference genome map
        posAllelePackMap = new HashMap<>(vlEndIndex-vlStartIndex);  // the map key is position, what is the value data type?
        positions = new int[vlEndIndex-vlStartIndex];   // ok, it is the same. but how?
        for (int i = vlStartIndex; i < vlEndIndex; i++) {
            posRefMap.put(vl.positions[i], String.valueOf(frb.getBase(vl.positions[i]-1))); // base in reference genome, but why the positionIndex needs to - 1? every base needs a K-V map, is big for memory.
            posAllelePackMap.put(vl.positions[i], vl.getAllelePacks(i));    // key: position; value: AllelePackage[] (alleles at one position)
            positions[i-vlStartIndex] = vl.positions[i];    // make positions in this region store into a new int[]
        }
        Set<String> taxaSet = taxaBamsMap.keySet(); // created in retrieveAppParameters. key: taxa; value: List<String>
        ArrayList<String> taxaList = new ArrayList(taxaSet);    // from set to list? why directly to list?
        Collections.sort(taxaList); // the super-super-class of ArrayList implements Collection
        Dyad<int[][], int[]> d = FastCall2.getBins(this.regionStart, this.regionEnd, FastCall2.scanBinSize);    // 5M bp bin parallel process tradition
        int[][] binBound = d.getFirstElement();
        int[] binStarts = d.getSecondElement();
        LongAdder counter = new LongAdder();
        ExecutorService pool = Executors.newFixedThreadPool(this.threadsNum);
        List<Future<IndiCount>> resultList = new ArrayList<>();
        for (int i = 0; i < taxaList.size(); i++) { // recursive for each taxon
            List<String> bamPaths = taxaBamsMap.get(taxaList.get(i));   // a taxon can have many bam, is it the same in disc?
            StringBuilder sb = new StringBuilder(samtoolsPath); // sb start with samtools path
            sb.append(" mpileup --no-output-ends ").append(this.baqMode).append("-q ").append(this.mappingQThresh).append(" -Q ").append(this.baseQThresh).append(" -f ").append(this.referenceFileS);
            for (int j = 0; j < bamPaths.size(); j++) {
                sb.append(" ").append(bamPaths.get(j));
            }
            sb.append(" -l ").append(vLibPosFileS).append(" -r ");
            sb.append(chrom).append(":").append(this.regionStart).append("-").append(this.regionEnd);
            String command = sb.toString();
//            System.out.println(command);
            IndiCount idv = new IndiCount(command, taxaList.get(i), binBound, binStarts, bamPaths, counter);    // return an indiCount class and write .iac files
            Future<IndiCount> result = pool.submit(idv);    // why there is a pool and need to do the for loop?
            resultList.add(result); // is there a need to do this?
        }
        try {
            pool.shutdown();
            pool.awaitTermination(Long.MAX_VALUE, TimeUnit.SECONDS);
        }
        catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }

    }

    class IndiCount implements Callable<IndiCount> {
        String command = null;
        String taxonName = null;
        String indiTaxonDirS = null;
        int[][] binBound = null;
        int[] binStarts = null;
        List<String> bamPaths = null;
        LongAdder counter = null;
        DataOutputStream dos = null;
        int currentBinIndex = Integer.MIN_VALUE;

        public IndiCount (String command, String taxonName, int[][] binBound, int[] binStarts, List<String> bamPaths, LongAdder counter) {
            this.command = command;
            this.taxonName = taxonName;
            this.binBound = binBound;
            this.binStarts = binStarts;
            this.bamPaths = bamPaths;
            this.counter = counter;
            String indiCountFolderS = new File(outputDirS, subDirS[1]).getAbsolutePath();
            indiTaxonDirS = new File(indiCountFolderS, taxonName).getAbsolutePath();
            new File (indiTaxonDirS).mkdir();
        }

        public void closeDos () {
            try {
                dos.flush();
                dos.close();
            }
            catch (Exception e) {
                e.printStackTrace();
                System.exit(1);
            }
        }

        public void setDos (int queryPos) {
            int binIndex = Arrays.binarySearch(binStarts, queryPos);    // find the bin
            if (binIndex < 0) binIndex = -binIndex-2;   // to the index of the matched position in this bin
            if (binIndex != currentBinIndex) {  // every bin write a
                if (currentBinIndex > -1) this.closeDos();
                StringBuilder sb = new StringBuilder();
                sb.append(chrom).append("_").append(binBound[binIndex][0]).append("_").append(binBound[binIndex][1]).append(".iac.gz");
                String outfileS = new File (indiTaxonDirS, sb.toString()).getAbsolutePath();
                dos = IOUtils.getBinaryGzipWriter(outfileS);
                try {
                    dos.writeUTF(this.taxonName);
                    dos.writeShort((short)chrom);
                    dos.writeInt(binBound[binIndex][0]);
                    dos.writeInt(binBound[binIndex][1]);
                    vlBinStartIndex = vl.getStartIndex(binBound[binIndex][0]);
                    vlBinEndIndex = vl.getEndIndex(binBound[binIndex][1]-1);    // output the position index that <= bin end
//                    if (vlBinStartIndex < 0 || vlBinEndIndex < 0) {
////                        unlikely to happen
//                    }
                    dos.writeInt(vlBinEndIndex-vlBinStartIndex);    // positions number
                }
                catch (Exception e) {
                    e.printStackTrace();
                    System.exit(1);
                }
                currentBinIndex = binIndex;
            }
        }

        public void writeAlleleCounts (int[] alleleCounts) {
            try {
                dos.writeByte((byte)alleleCounts.length);   // allele number
                for (int i = 0; i < alleleCounts.length; i++) {
                    if (alleleCounts[i] < Short.MAX_VALUE) {
                        dos.writeShort((short)alleleCounts[i]);
                    }
                    else {
                        dos.writeShort(Short.MAX_VALUE);
                    }
                }
            }
            catch (Exception e) {
                e.printStackTrace();
                System.exit(1);
            }
        }

        public void writeMissing () {
            try {
                dos.writeByte((byte)-1);    // means nothing? so what is something?
            }
            catch (Exception e) {
                e.printStackTrace();
                System.exit(1);
            }
        }

        public void writeEmptyFiles() {
            StringBuilder sb = new StringBuilder();
            for (int i = 0; i < binStarts.length; i++) {
                sb.setLength(0);
                sb.append(chrom).append("_").append(binBound[i][0]).append("_").append(binBound[i][1]).append(".iac.gz");
                String outfileS = new File (indiTaxonDirS, sb.toString()).getAbsolutePath();
                File outf = new File (outfileS);
                if (outf.exists()) continue;
                try {
                    dos = IOUtils.getBinaryGzipWriter(outfileS);
                    dos.writeUTF(this.taxonName);
                    dos.writeShort((short)chrom);
                    dos.writeInt(binBound[i][0]);
                    dos.writeInt(binBound[i][1]);
                    vlBinStartIndex = vl.getStartIndex(binBound[i][0]);
                    vlBinEndIndex = vl.getEndIndex(binBound[i][1]-1);
                    int len = vlBinEndIndex-vlBinStartIndex;
                    if (vlBinStartIndex == Integer.MIN_VALUE || vlBinEndIndex == Integer.MIN_VALUE) {
                        len = 0;
                    }
                    dos.writeInt(len);
                    for (int j = 0; j < len; j++) {
                        this.writeMissing();
                    }
                    this.closeDos();
                }
                catch (Exception e) {
                    e.printStackTrace();
                }
            }
        }

        @Override
        public IndiCount call() throws Exception {
            try {
                Runtime rt = Runtime.getRuntime();
                Process p = rt.exec(command);   // directly run?
                String temp = null;
                BufferedReader br = new BufferedReader(new InputStreamReader(p.getInputStream()));  // why input stream is from process?
                BufferedReader bre = new BufferedReader(new InputStreamReader(p.getErrorStream())); // what error? the error after launch?
                DataOutputStream dis = null;    // wrong written
                String current = br.readLine(); // what is the content in br?
                List<String> currentList = null;
                int currentPosition = -1;
                if (current != null) {
                    currentList = PStringUtils.fastSplit(current);  // the use of List<String>, can learn more about PStringUtils
                    currentPosition = Integer.parseInt(currentList.get(1)); // prove br is the output of mpileup
                }
                StringBuilder baseS = new StringBuilder();  // why base needs a StringBuilder?
                StringBuilder indelSb = new StringBuilder();
                for (int i = 0; i < positions.length; i++) {
                    this.setDos(positions[i]);  // write the basic information for each bin (identify by binIndex)
                    //at the end of pileup file
                    if (current == null) {
                        this.writeMissing();    // no variation for this taxon at this position
                    }
                    else {
                        if (positions[i] == currentPosition) {  // so one call only record one currentPosition? but is mpileup run for every taxon and for each bin? it reruns everything. there must be a difference.
                            String ref = posRefMap.get(currentPosition);    // get ref base, but no use?
                            AllelePackage[] altAlleles = posAllelePackMap.get(currentPosition); // AllelePackage[] from VariationLibrary
                            baseS.setLength(0);
                            int siteDepth = 0;
                            for (int j = 0; j < bamPaths.size(); j++) {
                                siteDepth+=Integer.parseInt(currentList.get(3+j*3));    // count whole depth; also for the same taxon; no difference with disc. but disc did not record the depth, it directly decides the heterozygous or homozygous
                                baseS.append(currentList.get(4+j*3));   // baseS means base String
                            }
                            if (siteDepth == 0) {
                                this.writeMissing();
                            }
                            else {
//                                FastCall2.removeFirstPositionSign(baseS);
                                int[] alleleCounts = getAlleleCounts (altAlleles, baseS.toString().toUpperCase(), siteDepth, indelSb);  // return an int[]. it is only the allele base count number (with reference genotype at index 0)
                                // why not directly store the count information at the first time when doing mpileup?
                                this.writeAlleleCounts(alleleCounts);   // same question: why not write at the first time? no built allelePack? but the count information can be stored in a file for future count?
                            }
                            current = br.readLine();
                            if (current != null) {
                                currentList = PStringUtils.fastSplit(current);
                                currentPosition = Integer.parseInt(currentList.get(1));
                            }
                        }
                        else if (positions[i] < currentPosition) {
                            this.writeMissing();
                        }
                        else {
                            System.out.println("Current position is greater than pileup position. It should not happen. Program quits");    // debug information
                            System.exit(1);
                        }
                    }
                }
                this.closeDos();
                br.close();
                p.waitFor();
                this.writeEmptyFiles();
                StringBuilder sb = new StringBuilder();
                sb.append("mpileup command: ").append(command).append("\n");    // record command in the end
                sb.append("mpileup error profile: ");
                while ((current = bre.readLine()) != null) {
                    sb.append(current).append("\n");
                }
                sb.append("Individual allele counting is completed for taxon ").append(this.taxonName).append("\n");
                System.out.println(sb.toString());

            }
            catch (Exception ee) {
                ee.printStackTrace();
            }
            counter.increment();
            int cnt = counter.intValue();   // the number of taxa counted
            if (cnt%50 == 0) System.out.println("Finished individual genotype allele counting in " + String.valueOf(cnt) + " taxa. Total: " + String.valueOf(taxaBamsMap.size()));
            return this;
        }
    }

    private void deleteTemperateFile () {
        File f1 = new File(outputDirS, subDirS[0]);
        File f2 = new File(outputDirS, subDirS[1]);
        try {
            FileUtils.cleanDirectory(new File(outputDirS, subDirS[0]));
            FileUtils.cleanDirectory(new File(outputDirS, subDirS[1]));
        }
        catch (Exception e) {
            e.printStackTrace();
        }
        f1.delete();
        f2.delete();
        new File (vLibPosFileS).delete();
    }

    private String getInfoAndGenotypes (List<short[]> siteCountList, AllelePackage[] altAlleles) {
        StringBuilder genoSB = new StringBuilder("GT:AD:GL");
        StringBuilder infoSB = new StringBuilder();
        int dp = 0; // DP: Total Depth
        int nz = 0; // NZ: Number of taxa with called genotypes
        int alleleNumber = altAlleles.length+1; // includes ref
        int[] adCnt = new int[alleleNumber];    // AD: Total allelelic depths in order listed starting with REF
        int[] acCnt = new int[alleleNumber];    // AC: Numbers of allele occurence across taxa in order listed
        int[][] gnCnt = new int[alleleNumber][alleleNumber];    // GN: Number of taxa with genotypes AA,AB,BB or AA,AB,AC,BB,BC,CC if 2 alt alleles
        int ht = 0; // HT: Number of heterozygotes
        for (int i = 0; i < siteCountList.size(); i++) {    // recur for each taxon, including ref
            if (siteCountList.get(i) == null) { // no variation, therefore not called. but how to decide the depth? the mean depth of alignment? **
                nz++;   // called genotype is that same as ref? is the number of taxa
                genoSB.append("\t./.");
                continue;
            }
            for (int j = 0; j < alleleNumber; j++) {
                int currentCount = siteCountList.get(i)[j]; // depth for taxon i and genotype j
                dp+=currentCount;   // total depth
                adCnt[j] += currentCount;
            }
            String genoS = this.getGenotypeByShort(siteCountList.get(i));   // count step for taxa i *
            genoSB.append("\t").append(genoS);
            int g1 = genoS.charAt(0)-'0';   // genotype 1
            int g2 = genoS.charAt(2)-'0';   // genotype 2
            acCnt[g1]++;
            acCnt[g2]++;
            gnCnt[g1][g2]++;
            if (g1 != g2) ht++;
        }
        nz = siteCountList.size() - nz; // real NZ?
        int sum = 0;
        for (int i = 0; i < acCnt.length; i++) {
            sum+=acCnt[i];  // sum all the allele count
        }
        float maf = (float)((double)acCnt[0]/sum);  // choose one and count the frequency
        if (maf>0.5) maf = (float)(1-maf);  // find the smaller one
        infoSB.append("DP=").append(dp).append(";NZ=").append(nz).append(";AD=");
        for (int i = 0; i < adCnt.length; i++) {    // for each genotype
            infoSB.append(adCnt[i]).append(",");
        }
        infoSB.deleteCharAt(infoSB.length()-1);
        infoSB.append(";AC=");  // sum of AC is 2*NZ
        for (int i = 0; i < acCnt.length; i++) {
            infoSB.append(acCnt[i]).append(",");
        }
        infoSB.deleteCharAt(infoSB.length()-1);
        infoSB.append(";IS=");
        for (int i = 0; i < altAlleles.length; i++) {
            infoSB.append(altAlleles[i].getIndelSeq()).append(",");
        }
        infoSB.deleteCharAt(infoSB.length()-1);
        infoSB.append(";GN=");  // genotype number, sum of GN is NZ
        for (int i = 0; i < gnCnt.length; i++) {
            for (int j = i; j < gnCnt.length; j++) {
                infoSB.append(gnCnt[i][j]).append(",");
            }
        }
        infoSB.deleteCharAt(infoSB.length()-1);
        infoSB.append(";HT=").append(ht).append(";MAF=").append(maf);   // output HT, which is the mid-number at GN
        infoSB.append("\t").append(genoSB);
        return infoSB.toString();
    }

    private int[] getAlleleCounts (AllelePackage[] altAlleles, String baseS, int siteDepth, StringBuilder indelSb) {
        byte[] baseB = baseS.getBytes();
        int[] altAlleleCounts = new int[altAlleles.length]; // the count is in one taxon?
        int index = Integer.MIN_VALUE;
        int vCnt = 0;
        for (int i = 0; i < baseB.length; i++) {
            byte queryAlleleCoding = FastCall2.pileupAscIIToAlleleCodingMap.get(baseB[i]);  // from byte (AscIIs) to byte ({000, 001, 010, 011, 100, 101})
            int queryIndelLength = 0;
            index = Arrays.binarySearch(AlleleEncoder.alleleCodings, queryAlleleCoding);    // find indel
            if (index > 3) {    // insertion and deletion are {4, 5}
                int startIndex = i+1;
                int endIndex = i+2;
                for (int j = i+2; j < baseB.length; j++) {
                    if (baseB[j] > 57) {
                        endIndex = j;
                        break;
                    }
                }
                indelSb.setLength(0);
                for (int j = startIndex; j < endIndex; j++) {
                    indelSb.append((char)baseB[j]);
                }
                queryIndelLength = Integer.parseInt(indelSb.toString());
                i+=indelSb.length();
                i+=queryIndelLength;
            }
            for (int j = 0; j < altAlleles.length; j++) {
                if (altAlleles[j].getAlleleCoding() == queryAlleleCoding && altAlleles[j].getIndelLength() == queryIndelLength) {
                    altAlleleCounts[j]++;   // count every base?
                    vCnt++;
                }
            }
        }
        int[] alleleCounts = new int[altAlleles.length+1];
        alleleCounts[0] = siteDepth - vCnt; // reference genotype count
        for (int i = 0; i < altAlleles.length; i++) {
            alleleCounts[i+1] = altAlleleCounts[i];
        }
        return alleleCounts;    // return a int[]
    }

    private String getGenotypeByShort (short[] cnt) {   // is for a single taxon
        //in case some allele depth is greater than maxFactorial, to keep the original allele counts
        short[] oriCnt = null;
        int n = cnt.length*(cnt.length+1)/2;    // for what?
        int[] likelihood = new int[n];
        int sum = 0;
        for (int i = 0; i < cnt.length; i++) sum+=cnt[i];   // sum of depths
        if (sum == 0) return "./."; // no mapped depth. by this do not need to plus the NZ?
        else if (sum > this.maxFactorial) { // the use of factorial
            oriCnt = new short[cnt.length]; // faster?
            System.arraycopy(cnt, 0, oriCnt, 0, cnt.length);    // copy all information to store original allele count information
            double portion = (double)this.maxFactorial/sum;
            for (int i = 0; i < cnt.length; i++) {
                cnt[i] = (short)(cnt[i]*portion);   // use portion as the number, change cnt
            }
            sum = this.maxFactorial;    // redefine the sum
        }
        // what is the algorithm??
        double coe = this.factorialMap.get(sum);    // use factorial for what? count what?
        for (int i = 0; i < cnt.length; i++) coe = coe/this.factorialMap.get(cnt[i]);   // what is the function? likelihood function? multiple probabilities together?
        double max = Double.MAX_VALUE;  // for what?
        int a1 = 0;
        int a2 = 0;
        for (int i = 0; i < cnt.length; i++) {  // every allele
            for (int j = i; j < cnt.length; j++) {  // another allele (including itself)
                int index = (j*(j+1)/2)+i;  // function?
                double value = Double.MAX_VALUE;
                if (i == j) {
                    value = -Math.log10(coe*Math.pow((1-0.75*this.combinedErrorRate), cnt[i])*Math.pow(this.combinedErrorRate /4, (sum-cnt[i])));
                    // what is the function? and the theory?
                    // to count the likelihood, so what is the distribution?
                }
                else {
                    value = -Math.log10(coe*Math.pow((0.5-this.combinedErrorRate /4), cnt[i]+cnt[j])*Math.pow(this.combinedErrorRate /4, (sum-cnt[i]-cnt[j])));
                }
                if (value < max) {  // get the smallest one (why the minimum?)
                    max = value;
                    a1 = i;
                    a2 = j;
                }
                likelihood[index] = (int)Math.round(value); // finish
                // actually, the memory can store all of these data
            }
        }
        StringBuilder sb = new StringBuilder();
        sb.append(a1).append("/").append(a2).append(":");   // GT
        // AD
        if (sum > this.maxFactorial) {
            for (int i = 0; i < oriCnt.length; i++) sb.append(oriCnt[i]).append(",");
        }
        else {
            for (int i = 0; i < cnt.length; i++) sb.append(cnt[i]).append(",");
        }
        sb.deleteCharAt(sb.length()-1); sb.append(":"); // separator
        for (int i = 0; i < likelihood.length; i++) sb.append(likelihood[i]).append(",");   // likelihood
        sb.deleteCharAt(sb.length()-1);
        return sb.toString();
    }

    private void processVariationLibrary () {
        StringBuilder sb = new StringBuilder();
        sb.append(this.chrom).append("_").append(this.regionStart).append("_").append(regionEnd).append(".pos.txt");    // new file
        this.vLibPosFileS = new File (this.outputDirS, sb.toString()).getAbsolutePath();
        this.vl = new VariationLibrary(this.libFileS);  // read the VariationLibrary
        if (this.chrom != vl.getChrom()) {
            System.out.println("The chromosome number of library and the specified one do not match. Program quits.");
            System.exit(0);
        }
        try {
            vlStartIndex = vl.getStartIndex(this.regionStart);  // return start position index
            vlEndIndex = vl.getEndIndex(this.regionEnd);        // return end position index
            if (vlStartIndex == Integer.MIN_VALUE || vlEndIndex == Integer.MIN_VALUE) {
                System.out.println("The chromosome region was incorrectly set. Program quits.");
                System.exit(0);
            }
            BufferedWriter bw = IOUtils.getTextWriter(this.vLibPosFileS);
            for (int i = vlStartIndex; i < vlEndIndex; i++) {
                sb.setLength(0);
                sb.append(this.chrom).append("\t").append(vl.getPosition(i));   // write to position file
                bw.write(sb.toString());
                bw.newLine();
            }
            bw.flush();
            bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }

    private void creatFactorialMap () {
        this.factorialMap = HashIntDoubleMaps.getDefaultFactory().newMutableMap();  // why factorial? is the factorial like x!
        for (int i = 0; i < this.maxFactorial+1; i++) { // from 0 to 150
            this.factorialMap.put(i, factorial(i)); // longFactorials = new long[]{1L, 1L, 2L, 6L, 24L, 120L...
            // why use factorial map?
        }
    }

    public void mkDir () {
        File f = new File (this.outputDirS);
        f.mkdir();
        for (int i = 0; i < subDirS.length; i++) {
            f = new File(outputDirS, subDirS[i]);   // create 3 sub-directions {"indiVCF", "indiCounts", "VCF"}
            f.mkdir();
        }
    }

    private void parseTaxaBamMap(String taxaBamMapFileS) {
        this.taxaBamsMap = new HashMap<>();
        this.taxaCoverageMap = new HashMap<>();
        try {
            BufferedReader br = IOUtils.getTextReader(taxaBamMapFileS); // read file
            String temp = br.readLine();    // read the first column name row?
            ArrayList<String> taxaList = new ArrayList();   // length changeable
            ArrayList<String> pathList = new ArrayList();   // no use after
            int nBam = 0;
            while ((temp = br.readLine()) != null) {
                // duplicate with the read process in DiscoverVariation
                String[] tem = temp.split("\t");    // is there a need to set the length of String[]?
                taxaList.add(tem[0]);
                String[] bams = new String[tem.length-2] ;  // minus the taxa and coverage column
                for (int i = 0; i < bams.length; i++) {
                    bams[i] = tem[i+2]; // add bam path
                    pathList.add(bams[i]);  // why all the bam (with different taxon) are stored in pathList?
                }
                Arrays.sort(bams);
                List<String> bamList = Arrays.asList(bams); // Array -> List, transfer from bams
                taxaBamsMap.put(tem[0], bamList);
                taxaCoverageMap.put(tem[0], Double.valueOf(tem[1]));    // no usage?
                nBam+=bams.length;  // why not directly use pathList's length?
            }
            taxaNames = taxaList.toArray(new String[taxaList.size()]);  // List -> Array
            Arrays.sort(taxaNames);
            HashSet<String> taxaSet = new HashSet<>(taxaList);  // no repeat
            if (taxaSet.size() != taxaNames.length) {
                System.out.println("Taxa names are not unique. Programs quits");
                System.exit(0);
            }
            System.out.println("Created TaxaBamMap from" + taxaBamMapFileS);
            System.out.println("Taxa number:\t"+String.valueOf(taxaNames.length));
            System.out.println("Bam file number in TaxaBamMap:\t"+String.valueOf(nBam));
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
}

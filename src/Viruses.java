import org.biojava.nbio.core.exceptions.CompoundNotFoundException;
import org.biojava.nbio.core.sequence.DNASequence;
import org.biojava.nbio.core.sequence.template.AbstractSequence;

import java.io.IOException;
import java.util.Map;
import java.util.HashMap;
import java.util.Arrays;
import java.util.ArrayList;
import java.util.List;
import java.util.Comparator;
import java.util.stream.Collectors;

/**
 * Used example:
 * https://github.com/biojava/biojava-tutorial/blob/master/core/README.md
 */

public class Viruses {

    private static int number = 0;
    private static double[] ATG = {0, 0, 0, 0, 0, 0, 0, 0};
    private static double[] TAA = {0, 0, 0, 0, 0, 0, 0, 0};
    private static double[] TAG = {0, 0, 0, 0, 0, 0, 0, 0};
    private static double[] TGA = {0, 0, 0, 0, 0, 0, 0, 0};
    private static Map<String, Integer>[] diCodonList = new HashMap[8];
    private static final int[] diCodonListSum = {0, 0, 0, 0, 0, 0, 0, 0};
    private static Map<String, Integer>[] sortedDiCodons = new HashMap[8];

    public static void main(String[] args) throws Exception {
        execute("bacterial1.fasta", number);
        number++;
        execute("bacterial2.fasta", number);
        number++;
        execute("bacterial3.fasta", number);
        number++;
        execute("bacterial4.fasta", number);
        number++;
        execute("mamalian1.fasta", number);
        number++;
        execute("mamalian2.fasta", number);
        number++;
        execute("mamalian3.fasta", number);
        number++;
        execute("mamalian4.fasta", number);
        number++;

//        System.out.println("ATG: " + Arrays.toString(ATG));
//        System.out.println("TAA: " + Arrays.toString(TAA));
//        System.out.println("TAG: " + Arrays.toString(TAG));
//        System.out.println("TGA: " + Arrays.toString(TGA));
//        System.out.println(Arrays.toString(diCodonList));
//        System.out.println(Arrays.toString(diCodonListSum));
//        System.out.println(Arrays.toString(sortedDiCodons));
//        System.out.println();

//        5 užduotis
        calculateCodonMatrix();
        calculateDiCodonMatrix();

        System.out.println("Program is closing...");
    }

    private static List<DNASequence> findCodonSequences(List<DNASequence> seq) throws CompoundNotFoundException {
        List<DNASequence> codonList = new ArrayList<>();
        int i = 0;

        while (i < seq.size()) {
            if (seq.get(i).getSequenceAsString().equals("ATG")) {
                final int startPos = i;
                int j = i;

                while (j < seq.size()) {
                    if (Arrays.asList("TAA", "TAG", "TGA").contains(seq.get(j).getSequenceAsString())) {

                        List<String> filteredSeq = seq.subList(startPos, j + 1)
                                .stream()
                                .map(AbstractSequence::getSequenceAsString)
                                .collect(Collectors.toList());

                        codonList.add(new DNASequence(String.join("", filteredSeq)));
                        i = j;

                        break;
                    }
                    j++;
                }
            }
            i++;
        }

        return codonList;
    }

    private static void printSequenceList(String text, List<DNASequence> seq) {
        System.out.println(text);

        seq.forEach(System.out::println);

        System.out.println();
    }

    private static DNASequence findLongestCodonSequence(List<DNASequence> seq) {
        return seq
                .stream()
                .max(Comparator.comparingInt(AbstractSequence::getLength))
                .orElse(new DNASequence());
    }

    private static List<DNASequence> removeUselessCodons(List<DNASequence> seq) {
        return seq
                .stream()
                .filter(element -> element.getLength() >= 100)
                .collect(Collectors.toList());
    }

    private static void execute(String fileName, int number) throws IOException, CompoundNotFoundException {
        DNASequence data = Utils.readFastaFile(fileName);

        System.out.println("Executing: " + fileName);

//        1 užduotis
        System.out.println("1st exercise");
        int split = 3;
        List<DNASequence> listFrame1 = Utils.splitStringToFrames(0, data.getSequenceAsString(), split);
        List<DNASequence> listFrame2 = Utils.splitStringToFrames(1, data.getSequenceAsString(), split);
        List<DNASequence> listFrame3 = Utils.splitStringToFrames(2, data.getSequenceAsString(), split);

        List<DNASequence> reverseListFrame1 = Utils.splitStringToFrames(0, data.getReverseComplement().getSequenceAsString(), split);
        List<DNASequence> reverseListFrame2 = Utils.splitStringToFrames(1, data.getReverseComplement().getSequenceAsString(), split);
        List<DNASequence> reverseListFrame3 = Utils.splitStringToFrames(2, data.getReverseComplement().getSequenceAsString(), split);

        List<DNASequence> codonSequences = findCodonSequences(listFrame1);
        codonSequences.addAll(findCodonSequences(listFrame2));
        codonSequences.addAll(findCodonSequences(listFrame3));
        codonSequences.addAll(findCodonSequences(reverseListFrame1));
        codonSequences.addAll(findCodonSequences(reverseListFrame2));
        codonSequences.addAll(findCodonSequences(reverseListFrame3));
        printSequenceList("All frames combined are as follows:", codonSequences);

//        2 užduotis
        System.out.println("2nd exercise");
        DNASequence longestCodon = findLongestCodonSequence(codonSequences);
        System.out.println(longestCodon);

//        3 užduotis
        System.out.println("3rd exercise");
        List<DNASequence> optimizedCodons = removeUselessCodons(codonSequences);
        printSequenceList("Longest codons:", optimizedCodons);

//        4 užduotis
        System.out.println("4th exercise");
        findCodonSequence(optimizedCodons);

//        System.out.printf("Normalized ATG: %s%n", ATG[number]);
//        System.out.printf("Normalized TAA: %s%n", TAA[number]);
//        System.out.printf("Normalized TAG: %s%n", TAG[number]);
//        System.out.printf("Normalized TGA: %s%n", TGA[number]);

        findDiCodonSequence(optimizedCodons);
    }

    private static void findCodonSequence(List<DNASequence> seq) {
        String splitRegex = "(?<=\\G...)";

        String oneString = String.join("", Utils.dnaSequenceToString(seq));
        System.out.println(oneString);

        List<String> splitSeq = Arrays.asList(oneString.split(splitRegex));
//        splitSeq.forEach(System.out::println);

        Map<String, Integer> occurrences = Utils.countOccurrences(splitSeq);
        System.out.println("Codon sequences");
        occurrences.forEach((key, value) -> System.out.println("Sequence: " + key + ", frequency: " + value.toString()));

        int sum = occurrences.values()
                .stream()
                .mapToInt(Integer::intValue)
                .sum();
        System.out.println(sum);

        ATG[number] = occurrences.get("ATG") * 1.0 / sum;
        TAA[number] = occurrences.get("TAA") * 1.0 / sum;
        TAG[number] = occurrences.get("TAG") * 1.0 / sum;
        TGA[number] = occurrences.get("TGA") * 1.0 / sum;
    }

    private static void findDiCodonSequence(List<DNASequence> seq) throws CompoundNotFoundException {
        String text = String.join("", Utils.dnaSequenceToString(seq));

        DNASequence dna = new DNASequence(text);

        char[] translatedDna = Utils.getTranscriptionEngine()
                .translate(dna)
                .getSequenceAsString()
                .toCharArray();

        List<String> di = new ArrayList<>();

        int i = 0;
        while (i < translatedDna.length - 1) {
            if (translatedDna[i] != '*' &&
                    translatedDna[i + 1] != '*') {
                String stringBuilder =
                        String.valueOf(translatedDna[i]) + translatedDna[i + 1];
                di.add(stringBuilder);
            }
            i++;
        }

        Map<String, Integer> occurrences = Utils.countOccurrences(di);

        System.out.println("DiCodon sequences");
        occurrences.forEach((key, value) -> System.out.println("Sequence: " + key + ", frequency: " + value.toString()));

        diCodonList[number] = occurrences;

        diCodonListSum[number] = occurrences.values()
                .stream().mapToInt(Integer::intValue).sum();

        sortedDiCodons[number] = Utils.sortByValue(diCodonList[number]);
    }

    private static void calculateDiCodonMatrix() {
        int z1 = 0;
        double[][] distance = new double[number][number];

        while (z1 < number) {
            int z2 = 0;
            double[] distanceRow = new double[number];
            while (z2 < number) {
                double diCodonSum = 0;
                int z3 = 0;

                if (sortedDiCodons[z1].size() < sortedDiCodons[z2].size()) {
                    while (z3 < sortedDiCodons[z1].size()) {
                        diCodonSum += Math.pow(
                                getSpecificValueFromDictionary(z1, z3) * 1.0f / diCodonListSum[z1] -
                                        getSpecificValueFromDictionary(z2, z3) * 1.0f / diCodonListSum[z2], 2);
                        z3++;
                    }
                } else {
                    while (z3 < sortedDiCodons[z2].size()) {
                        diCodonSum += Math.pow(
                                getSpecificValueFromDictionary(z1, z3) * 1.0f / diCodonListSum[z1] -
                                        getSpecificValueFromDictionary(z2, z3) * 1.0f / diCodonListSum[z2], 2);
                        z3++;
                    }
                }

                distanceRow[z2] = diCodonSum;
                z2++;
            }
            distance[z1] = distanceRow;
            z1++;
        }

        System.out.println("DiCodon Matrix output");
        System.out.println("B1 | B2 | B3 | B4 | M1 | M2 | M3 | M4");
        for (double[] element : distance) {
            for (double value : element) {
                System.out.print(String.format("%.6f", value) + " ");
            }
            System.out.println();
        }
    }

    private static void calculateCodonMatrix() {
        int z1 = 0;
        double[][] distance = new double[number][number];

        while (z1 < number) {
            int z2 = 0;
            double[] distanceRow = new double[number];
            while (z2 < number) {
                distanceRow[z2] = Math.pow(ATG[z1] - ATG[z2], 2) +
                        Math.pow(TAA[z1] - TAA[z2], 2) +
                        Math.pow(TAG[z1] - TAG[z2], 2) +
                        Math.pow(TGA[z1] - TGA[z2], 2);
                z2++;
            }
            distance[z1] = distanceRow;
            z1++;
        }

        System.out.println("Codon Matrix output");
        System.out.println("B1 | B2 | B3 | B4 | M1 | M2 | M3 | M4");
        for (double[] element : distance) {
            for (double value : element) {
                System.out.print(String.format("%.6f", value) + " ");
            }
            System.out.println();
        }
    }

    private static int getSpecificValueFromDictionary(int z1, int z3) {
        List<Map.Entry<String, Integer>> collect = sortedDiCodons[z1].entrySet()
                .stream()
                .sorted(Map.Entry.comparingByKey())
                .collect(Collectors.toList());

        return collect.get(z3).getValue();
    }

}

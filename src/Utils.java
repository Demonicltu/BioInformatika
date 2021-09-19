import org.biojava.nbio.core.exceptions.CompoundNotFoundException;
import org.biojava.nbio.core.sequence.DNASequence;
import org.biojava.nbio.core.sequence.compound.AmbiguityDNACompoundSet;
import org.biojava.nbio.core.sequence.compound.AmbiguityRNACompoundSet;
import org.biojava.nbio.core.sequence.compound.NucleotideCompound;
import org.biojava.nbio.core.sequence.io.FastaReaderHelper;
import org.biojava.nbio.core.sequence.template.AbstractSequence;
import org.biojava.nbio.core.sequence.template.CompoundSet;
import org.biojava.nbio.core.sequence.transcription.TranscriptionEngine;

import java.io.File;
import java.io.IOException;
import java.util.Map;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.ArrayList;
import java.util.List;
import java.util.stream.Collectors;

public class Utils {

    public static DNASequence readFastaFile(String fileName) throws IOException {
        LinkedHashMap<String, DNASequence> dnaSequence =
                FastaReaderHelper.readFastaDNASequence(new File("inputs/" + fileName));

        return dnaSequence.values().iterator().next();
    }

    public static List<DNASequence> splitStringToFrames(int position, String data, int split) throws CompoundNotFoundException {
        StringBuilder stringBuilder = new StringBuilder();
        List<DNASequence> splitStringsList = new ArrayList<>();
        char[] chars = data.toCharArray();

        int pos = 0;
        for (int i = position; i < chars.length; i++) {
            if (pos >= 3) {
                pos = 0;
                splitStringsList.add(new DNASequence(stringBuilder.toString()));
                stringBuilder = new StringBuilder();
            }

            stringBuilder.append(chars[i]);
            pos++;
        }

        return splitStringsList;
    }

    public static Map<String, Integer> countOccurrences(List<String> splitSeq) {
        Map<String, Integer> occurences = new HashMap<>();

        splitSeq.stream()
                .collect(Collectors.groupingBy(s -> s))
                .forEach((k, v) -> occurences.put(k, v.size()));

        return occurences;
    }

    public static <K, V extends Comparable<? super V>> Map<K, V> sortByValue(Map<K, V> map) {
        List<Map.Entry<K, V>> list = new ArrayList<>(map.entrySet());
        list.sort(Map.Entry.comparingByValue());

        Map<K, V> result = new LinkedHashMap<>();
        for (Map.Entry<K, V> entry : list) {
            result.put(entry.getKey(), entry.getValue());
        }

        return result;
    }

    public static List<String> dnaSequenceToString(List<DNASequence> seq) {
        return seq.stream()
                .map(AbstractSequence::getSequenceAsString)
                .collect(Collectors.toList());
    }

    public static TranscriptionEngine getTranscriptionEngine() {
        AmbiguityDNACompoundSet ambiguityDNACompoundSet = AmbiguityDNACompoundSet.getDNACompoundSet();
        CompoundSet<NucleotideCompound> nucleotideCompoundSet = AmbiguityRNACompoundSet.getRNACompoundSet();

        return new TranscriptionEngine.Builder()
                .dnaCompounds(ambiguityDNACompoundSet)
                .rnaCompounds(nucleotideCompoundSet)
                .build();
    }

}

package patric;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;

/**
 *
 * @author tobias
 */
public class Trna {

    private final HashMap<String, HashSet<String>> code_table = new HashMap();
    private final HashMap<String, HashSet<String>> aa_tRNA = new HashMap();
    private final HashMap<Character, Character> codon_pairs = new HashMap();
    private final HashMap<String, Float> codon_usage = new HashMap();
    private final HashSet<String> start_codons = new HashSet(Arrays.asList(new String[]{"UUG","CUG","AUU","AUC","AUA","AUG","GUG"}));

    public Trna() {
        initCodons();
    }

    public Trna(String path, Strain strain) throws IOException {
        initCodons();
        readCodonUsage(path + "codon_usage_", strain);
        readCodeTable(path + "genetic_code_table.txt");
        read_tRNA(path + "tRNA_" + strain + ".txt");
        toFile(path + "missing_aa_codons_" + strain.toString() + ".csv", strain);
    }

    private void initCodons() {
        codon_pairs.put('A', 'U');
        codon_pairs.put('T', 'A');
        codon_pairs.put('G', 'C');
        codon_pairs.put('C', 'G');
    }

    private void readCodonUsage(String path, Strain strain) throws IOException {
        BufferedReader br = new BufferedReader(new FileReader(path + strain.toString() + ".txt"));
        String line;
        while ((line = br.readLine()) != null) {
            String[] split = line.split("\t");
            for (String usage : split) {
                String[] codon_use = usage.split(" ");
                codon_usage.put(codon_use[0], Float.parseFloat(codon_use[1]));
            }
        }
        br.close();
    }

    private void readCodeTable(String file) throws IOException {
        BufferedReader br = new BufferedReader(new FileReader(file));
        String line;
        br.readLine();
        while ((line = br.readLine()) != null && !line.equals("")) {
            String[] split = line.split("\t");
            HashSet<String> codons = new HashSet<>();
            codons.addAll(Arrays.asList(split[2].split(" ")));
            code_table.put(split[1], codons);
        }
        br.close();
    }

    private void read_tRNA(String file) throws IOException {
        BufferedReader br = new BufferedReader(new FileReader(file));
        String line;
        br.readLine();
        while ((line = br.readLine()) != null && !line.equals("")) {
            String[] split = line.split("\t")[13].split("-");
            HashSet<String> set = aa_tRNA.containsKey(split[1]) ? aa_tRNA.get(split[1]) : new HashSet<String>();
            set.add(antiSense(split[2]));
            aa_tRNA.put(split[1], set);
        }
        br.close();
    }

    private String antiSense(String codon) {
        return new StringBuilder().append(codon_pairs.get(codon.charAt(2))).append(codon_pairs.get(codon.charAt(1))).append(codon_pairs.get(codon.charAt(0))).toString();
    }

    private void toFile(String path, Strain strain) throws FileNotFoundException {
        StringBuilder sb = new StringBuilder("Strain: " + strain.toString() + "\nMissing amino acids:\n");
        for (Map.Entry<String, HashSet<String>> entrySet : code_table.entrySet()) {
            String key = entrySet.getKey();
            if (!aa_tRNA.containsKey(key) && !key.equalsIgnoreCase("Stop")) {
                sb.append(key).append(", ");
            }
        }
        sb.append("\ncodons for existing amino acids\tmissing(usage)\tpresent(usage):\t ratio\n");
        for (Map.Entry<String, HashSet<String>> entrySet : code_table.entrySet()) {
            String key = entrySet.getKey();
            HashSet<String> v = entrySet.getValue();
            if (aa_tRNA.containsKey(key)) {
                HashSet<String> missing = new HashSet<>(v);
                HashSet<String> present = new HashSet<>(v);
                missing.removeAll(aa_tRNA.get(key));
                present.retainAll(aa_tRNA.get(key));
                if (v.size() > 0) {
                    sb.append(key).append('\t');
                    float sum_missing = 0, sum_present = 0;
                    for (String missing_codon : missing) {
                        String codon = start_codons.contains(missing_codon)? missing_codon : missing_codon.toLowerCase();
                        sb.append(codon).append('(').append(codon_usage.get(missing_codon)).append("), ");
                        sum_missing += codon_usage.get(missing_codon);
                    }
                    sb.append('\t');
                    for (String present_codon : present) {
                        String codon = start_codons.contains(present_codon)? present_codon : present_codon.toLowerCase();
                        sb.append(codon).append('(').append(codon_usage.get(present_codon)).append("), ");
                        sum_present += codon_usage.get(present_codon);
                    }
                    sb.append('\t').append(sum_missing/sum_present).append('\n');
                }
            }
        }
        //Printwriter, to file ...
        PrintWriter pw = new PrintWriter(path);
        pw.write(sb.toString());
        pw.close();
    }

    public static void main(String[] args) throws IOException {
        String path = "/home/h/harrert/Dropbox/UNI/BACHELOR/Daten_Ergebnisse/tRNA/"; //laptop: "/home/tobias/Coxiella/tRNA/"
        Trna t_177 = new Trna(path, Strain.Q177);
        Trna t_154 = new Trna(path, Strain.Q154);
        Trna t_493 = new Trna(path, Strain.rsa493);
    }
}

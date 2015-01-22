package patric;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;

/**
 *
 * @author tobias
 */
public class Mapper {

    private final HashMap<String, String> locus_group = new HashMap();
    private final HashMap<String, boolean[]> group_genome = new HashMap();
    private final HashMap<String, String> figFam_locus_493 = new HashMap();//82552
    private final HashMap<String, String> figFam_locus_177 = new HashMap();//107188
    private final HashMap<String, String> figFam_locus_154 = new HashMap();//77120
    private final HashMap<String, String> figFam_function = new HashMap();
    private final HashMap<String, String> patric_refseq = new HashMap();

    public Mapper(String path, String outpath,boolean ommit) throws IOException {
        map_patric_refseq(path+"locus_patric_refseq_rsa493.txt");
        map_patric_refseq(path+"locus_patric_refseq_q177.txt");
        map_patric_refseq(path+"locus_patric_refseq_q154.txt");
        process_FigFam(path + "ProteinFamilyFeatures.txt");
        toFile(outpath, ommit);
    }
    
    private void map_patric_refseq(String file) throws IOException{
        BufferedReader br = new BufferedReader(new FileReader(file));
        String line;
        while ((line = br.readLine()) != null && !"".equals(line)){
            String[] split = line.split("\t");
            patric_refseq.put(split[2], split[3]);
        }
    }

    private void process_FigFam(String file) throws IOException {
        BufferedReader br = new BufferedReader(new FileReader(file));
        String line;
        br.readLine();
        while ((line = br.readLine()) != null && !line.equals("")) {
            String[] split = line.split("\t");//0: Group, 1:Genome , 3:Locus
            locus_group.put(patric_refseq.getOrDefault(split[3], split[3]), split[0]);
            figFam_function.putIfAbsent(split[0], split[10]);
            if (split[3].contains("82552")) {//RSA_493
                if (figFam_locus_493.containsKey(split[0])) {
                    figFam_locus_493.put(split[0], figFam_locus_493.get(split[0]) + "," + patric_refseq.getOrDefault(split[3], split[3]));//.split("_")[1] entfernt redundanten Patric_Id header: 82552_0123 -> 0123
                } else {
                    figFam_locus_493.put(split[0], patric_refseq.getOrDefault(split[3], split[3]));
                }
            } else if (split[3].contains("107188")) {//Q177
                if (figFam_locus_177.containsKey(split[0])) {
                    figFam_locus_177.put(split[0], figFam_locus_177.get(split[0]) + "," + patric_refseq.getOrDefault(split[3], split[3]));
                } else {
                    figFam_locus_177.put(split[0], patric_refseq.getOrDefault(split[3], split[3]));
                }
            } else if (split[3].contains("77120")) {//Q154
                if (figFam_locus_154.containsKey(split[0])) {
                    figFam_locus_154.put(split[0], figFam_locus_154.get(split[0]) + "," + patric_refseq.getOrDefault(split[3], split[3]));
                } else {
                    figFam_locus_154.put(split[0], patric_refseq.getOrDefault(split[3], split[3]));
                }
            }
            boolean[] genomes = (group_genome.containsKey(split[0])) ? group_genome.get(split[0]) : new boolean[3];//bit vector genomes: 0: Q177, 1: Q154, 2: RSA 493
            if (split[1].equalsIgnoreCase("Coxiella burnetii 'MSU Goat Q177'")) {
                genomes[0] = true;
            } else if (split[1].equalsIgnoreCase("Coxiella burnetii CbuK_Q154")) {
                genomes[1] = true;
            } else if (split[1].equalsIgnoreCase("Coxiella burnetii RSA 493")) {
                genomes[2] = true;
            }
            group_genome.put(split[0], genomes);
        }
        br.close();
    }

    private void toFile(String path, boolean ommit) throws FileNotFoundException {
        StringBuilder sb = new StringBuilder("FigFam\tQ177|Q154|RSA_493\tPatric Locus tags Q177 (107188_...)\tQ154 (77120_...)\tRSA_493 (82552_...)\tnumber of proteins\tortholog\tfunction\n");
        HashSet<String> groups = new HashSet();
        for (Map.Entry<String, String> entrySet : locus_group.entrySet()) {
            String group = entrySet.getValue();
            boolean[] genomes = group_genome.get(group);
            if (!groups.contains(group)) {
                int l_154 = figFam_locus_154.containsKey(group) ? figFam_locus_154.get(group).split(",").length : 0;
                int l_177 = figFam_locus_177.containsKey(group) ? figFam_locus_177.get(group).split(",").length : 0;
                int l_493 = figFam_locus_493.containsKey(group) ? figFam_locus_493.get(group).split(",").length : 0;
                int length = l_154 + l_177 + l_493;
                String function = figFam_function.get(group).contains("ypothetical protein") ? "---" : figFam_function.get(group);
                if (genomes[0] && genomes[1] && genomes[2]) {
                    int orthologs = Math.min(Math.min(l_154, l_177), l_493);
                    if (ommit) {
                        if (l_177 != l_154 || l_154 != l_493) {//length % orthologs != 0
                            sb.append(group).append('\t').append(genomes[0] ? "1" : "0").append(genomes[1] ? "1" : "0").append(genomes[2] ? "1" : "0").append('\t').append(figFam_locus_177.getOrDefault(group, "---")).append('\t').append(figFam_locus_154.getOrDefault(group, "---")).append('\t').append(figFam_locus_493.getOrDefault(group, "---")).append('\t').append(length).append(" (").append(l_177).append("|").append(l_154).append("|").append(l_493).append(")").append('\t').append(orthologs).append('\t').append(function).append('\n');
                        }
                    }
                    else{
                        sb.append(group).append('\t').append(genomes[0] ? "1" : "0").append(genomes[1] ? "1" : "0").append(genomes[2] ? "1" : "0").append('\t').append(figFam_locus_177.getOrDefault(group, "---")).append('\t').append(figFam_locus_154.getOrDefault(group, "---")).append('\t').append(figFam_locus_493.getOrDefault(group, "---")).append('\t').append(length).append(" (").append(l_177).append("|").append(l_154).append("|").append(l_493).append(")").append('\t').append(orthologs).append('\t').append(function).append('\n');
                    }
                } else if (!genomes[0] && genomes[1] && genomes[2]) {
                        sb.append(group).append('\t').append(genomes[0] ? "1" : "0").append(genomes[1] ? "1" : "0").append(genomes[2] ? "1" : "0").append('\t').append(figFam_locus_177.getOrDefault(group, "---")).append('\t').append(figFam_locus_154.getOrDefault(group, "---")).append('\t').append(figFam_locus_493.getOrDefault(group, "---")).append('\t').append(length).append(" (").append(l_177).append("|").append(l_154).append("|").append(l_493).append(")").append('\t').append(Math.min(l_154, l_493)).append('\t').append(function).append('\n');
                } else if (genomes[0] && !genomes[1] && genomes[2]) {
                        sb.append(group).append('\t').append(genomes[0] ? "1" : "0").append(genomes[1] ? "1" : "0").append(genomes[2] ? "1" : "0").append('\t').append(figFam_locus_177.getOrDefault(group, "---")).append('\t').append(figFam_locus_154.getOrDefault(group, "---")).append('\t').append(figFam_locus_493.getOrDefault(group, "---")).append('\t').append(length).append(" (").append(l_177).append("|").append(l_154).append("|").append(l_493).append(")").append('\t').append(Math.min(l_177, l_493)).append('\t').append(function).append('\n');
                } else if (genomes[0] && genomes[1] && !genomes[2]) {
                        sb.append(group).append('\t').append(genomes[0] ? "1" : "0").append(genomes[1] ? "1" : "0").append(genomes[2] ? "1" : "0").append('\t').append(figFam_locus_177.getOrDefault(group, "---")).append('\t').append(figFam_locus_154.getOrDefault(group, "---")).append('\t').append(figFam_locus_493.getOrDefault(group, "---")).append('\t').append(length).append(" (").append(l_177).append("|").append(l_154).append("|").append(l_493).append(")").append('\t').append(Math.min(l_154, l_177)).append('\t').append(function).append('\n');
                } else {
                    sb.append(group).append('\t').append(genomes[0] ? "1" : "0").append(genomes[1] ? "1" : "0").append(genomes[2] ? "1" : "0").append('\t').append(figFam_locus_177.getOrDefault(group, "---")).append('\t').append(figFam_locus_154.getOrDefault(group, "---")).append('\t').append(figFam_locus_493.getOrDefault(group, "---")).append('\t').append(length).append(" (").append(l_177).append("|").append(l_154).append("|").append(l_493).append(")").append("\t\t").append(function).append('\n');
                }
            }
            groups.add(group);
        }
        PrintWriter pw = new PrintWriter(path);
        pw.write(sb.toString());
        pw.close();
    }

    public static void main(String[] args) throws IOException {
        //Mapper m = new Mapper("/home/h/harrert/Coxiella/","/home/h/harrert/Desktop/figfam_genomes.csv");
        new Mapper("/home/tobias/Dropbox/UNI/BACHELOR/Daten_Ergebnisse/", "/home/tobias/Desktop/figfam_genomes.csv", false);
    }
}

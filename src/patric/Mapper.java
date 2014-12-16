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
    private final HashMap<String, String> patric_refSeq_locus = new HashMap();
    private final HashMap<String, String> figFam_refSeq_locus = new HashMap();

    public Mapper(String path) throws IOException{
        read_LocusTag("/home/tobias/Desktop/locus_patric_refseq_q177.txt");
        read_LocusTag("/home/tobias/Desktop/locus_patric_refseq_q154.txt");
        read_LocusTag("/home/tobias/Desktop/locus_patric_refseq_rsa493.txt");
        process_FigFam("/home/tobias/Desktop/ProteinFamilyFeatures.txt");
        toFile(path);
    }
    
    private void read_LocusTag(String file) throws IOException{
        BufferedReader br = new BufferedReader(new FileReader(file));
        String line;
        br.readLine();
        while ((line = br.readLine()) != null && !line.equals("")) {
            String[] split = line.split("\t");
            patric_refSeq_locus.put(split[2], split[3]);
        }
        br.close();
    }
    
    private void process_FigFam(String file) throws IOException {
        BufferedReader br = new BufferedReader(new FileReader(file));
        String line;
        br.readLine();
        while ((line = br.readLine()) != null && !line.equals("")) {
            String[] split = line.split("\t");//0: Group, 1:Genome , 3:Locus
            locus_group.put(split[3], split[0]);
            if(patric_refSeq_locus.containsKey(split[3])){
                if(figFam_refSeq_locus.containsKey(split[0])){
                    figFam_refSeq_locus.put(split[0], figFam_refSeq_locus.get(split[0])+","+patric_refSeq_locus.get(split[3]));
                }
                else{
                    figFam_refSeq_locus.put(split[0], patric_refSeq_locus.get(split[3]));
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
    
    private void toFile(String path) throws FileNotFoundException{
        StringBuilder sb = new StringBuilder("FigFam\tQ177\tQ154\tRSA_493\tRefSeq Locus tags\n");
        HashSet<String> groups = new HashSet();
        for (Map.Entry<String, String> entrySet : locus_group.entrySet()) {
            String locus = entrySet.getKey();
            String group = entrySet.getValue();
            boolean[] genomes = group_genome.get(group);
            if(!groups.contains(group)){
                sb.append(group).append('\t').append(genomes[0]? "1":"0").append('\t').append(genomes[1]?"1":"0").append('\t').append(genomes[2]?"1":"0").append('\t').append(figFam_refSeq_locus.getOrDefault(group,"no RefSeq locus mapped")).append('\n');
                
            }
            groups.add(group);
        }
        PrintWriter pw = new PrintWriter(path);
        pw.write(sb.toString());
        pw.close();
    }

    public static void main(String[] args) throws IOException {
        Mapper m = new Mapper("/home/tobias/Desktop/figfam_genomes.txt");
    }
}

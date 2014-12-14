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
    private final HashMap<String, Integer> group_families = new HashMap();

    private void processFile(String file) throws IOException {
        BufferedReader br = new BufferedReader(new FileReader(file));
        String line;
        br.readLine();
        while ((line = br.readLine()) != null && !line.equals("")) {
            String[] split = line.split("\t");//0: Group, 1:Genome , 3:Locus
            locus_group.put(split[3], split[0]);
            group_families.put(split[0], group_families.containsKey(split[0])? group_families.get(split[0])+1 : 1);
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
        StringBuilder sb = new StringBuilder("Locus Tag/FigFam\tQ177\tQ154\tRSA_493\n");
        HashSet<String> groups = new HashSet();
        for (Map.Entry<String, String> entrySet : locus_group.entrySet()) {
            String locus = entrySet.getKey();
            String group = entrySet.getValue();
            boolean[] genomes = group_genome.get(group);
            if(!groups.contains(group)){
                sb.append(group_families.get(group) == 1? locus : group).append('\t').append(genomes[0]).append('\t').append(genomes[1]).append('\t').append(genomes[2]).append('\n');
            }
            groups.add(group);
        }
        PrintWriter pw = new PrintWriter(path);
        pw.write(sb.toString());
        pw.close();
    }

    public static void main(String[] args) throws IOException {
        Mapper m = new Mapper();
        m.processFile("/home/tobias/Desktop/ProteinFamilyFeatures.txt");
        m.toFile("/home/tobias/Desktop/genomes_figfam_locus.txt");
    }
}

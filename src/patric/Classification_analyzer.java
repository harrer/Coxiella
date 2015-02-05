package patric;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import translation.Diff;

/**
 *
 * @author tobias
 */
public class Classification_analyzer {
    
    private static final HashMap<String, String> proteome = new HashMap<>();
    
    private static void read_proteome(String path) throws IOException{
        for (Strain strain : Strain.values()) {
            BufferedReader br = new BufferedReader(new FileReader(path+strain+".PATRIC.faa"));
            String line, header="";
            StringBuilder sb = new StringBuilder();
            while ((line = br.readLine()) !=null) {                
                if(line.startsWith(">")){
                    sb = new StringBuilder();
                    header = line.split("\\|")[3];
                }
                else if(!line.equals("")){
                    sb.append(line.trim());
                }
                else{
                    proteome.put(header, sb.toString());
                }
            }
            br.close();
        }
    }
    
    private static FIGFam[] process_FigFam(String file) throws IOException {
        BufferedReader br = new BufferedReader(new FileReader(file));
        HashMap<String, FIGFam> fams = new HashMap<>();
        String line;
        br.readLine();
        while ((line = br.readLine()) != null && !line.equals("")) {
            String[] split = line.split("\t");//0: Group, 1:Genome , 3:Locus
            FIGFam fam = fams.containsKey(split[0]) ? fams.get(split[0]) : new FIGFam(split[0]);
            if (split[3].contains("82552")) {//RSA_493
                fam.addMember(split[3], Strain.rsa493);
            } else if (split[3].contains("107188")) {//Q177
                if (Integer.parseInt(split[3].split("_")[1]) < 2159) {
                    fam.addMember(split[3], Strain.Q177);
                }
            } else if (split[3].contains("77120")) {//Q154
                fam.addMember(split[3], Strain.Q154);
            }
            fams.put(split[0], fam);
            /* boolean[] genomes = (group_genome.containsKey(split[0])) ? group_genome.get(split[0]) : new boolean[3];//bit vector genomes: 0: Q177, 1: Q154, 2: RSA 493
            if (split[1].equalsIgnoreCase("Coxiella burnetii 'MSU Goat Q177'")) {
            genomes[0] = true;
            } else if (split[1].equalsIgnoreCase("Coxiella burnetii CbuK_Q154")) {
            genomes[1] = true;
            } else if (split[1].equalsIgnoreCase("Coxiella burnetii RSA 493")) {
            genomes[2] = true;
            }
            group_genome.put(split[0], genomes);*/
        }
        br.close();
        return fams.values().toArray(new FIGFam[]{});
    }
    
    private static FIGFam[] find_specific_families(FIGFam[] figfams, Strain[] strains){
        ArrayList<FIGFam> fams = new ArrayList<>();
        for (FIGFam fam : figfams) {
            if(fam.hasStrain(strains, true)){
                fams.add(fam);
            }
        }
        return fams.toArray(new FIGFam[]{});
    }
    
    private static FIGFam[] find_orthologues(FIGFam[] figfams, Strain[] strains, double ratioThreshold, boolean equals){
        ArrayList<FIGFam> eqFams = new ArrayList<>();
        ArrayList<FIGFam> uneqFams = new ArrayList<>();
        int equal = 0, total = 0;
        for (FIGFam fam : figfams) {
            total++;
            if(fam.hasExactly_n_ProteinsPerMember(1)){
                //total++;
                String p1 = fam.getMembers().get(strains[0]).toArray(new String[]{})[0];
                String p2 = fam.getMembers().get(strains[1]).toArray(new String[]{})[0];
                if(Diff.matchingRatio(proteome.get(p1), proteome.get(p2)) >= ratioThreshold){
                    equal++;
                    eqFams.add(fam);
                }
                else{
                    uneqFams.add(fam);
                }
            }
            else{
                uneqFams.add(fam);
            }
        }
        System.out.println("eq:"+equal+"\ntotal: "+total);
        return equals? eqFams.toArray(new FIGFam[]{}) : uneqFams.toArray(new FIGFam[]{});
    }
    
    public static void main(String[] args) throws IOException {
        String path = "/home/tobias/Dropbox/UNI/BACHELOR/Daten_Ergebnisse/";
        read_proteome(path+"proteome/");
        FIGFam[] allFams = process_FigFam(path+"ProteinFamilyFeatures.txt");
        //Strain[] strains = new Strain[]{Strain.Q177, Strain.Q154};
        Strain[] strains = new Strain[]{Strain.Q154, Strain.rsa493};
        FIGFam[] specificFams = find_specific_families(allFams, strains);
        FIGFam[] equalFams = find_orthologues(specificFams, strains, 1, false);//false means unequal
        for (FIGFam equalFam : equalFams) {
            System.out.println(equalFam.getId());
        }
    }
}

package patric;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.function.Consumer;
import translation.Diff;
import lcs.LCS;

/**
 *
 * @author tobias
 */
public class Classification_analyzer {

    private static final HashMap<String, String> proteome = new HashMap<>();
    private static final HashMap<String, ArrayList<String>> seqMap = new HashMap<>();

    private static void read_proteome(String path) throws IOException {
        for (Strain strain : Strain.values()) {
            BufferedReader br = new BufferedReader(new FileReader(path + strain + ".PATRIC.faa"));
            String line, header = "";
            StringBuilder sb = new StringBuilder();
            while ((line = br.readLine()) != null) {
                if (line.startsWith(">")) {
                    sb = new StringBuilder();
                    header = line.split("\\|")[3];
                } else if (!line.equals("")) {
                    sb.append(line.trim());
                } else {
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
            FIGFam fam = fams.containsKey(split[0]) ? fams.get(split[0]) : new FIGFam(split[0], split[10]);
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

    private static FIGFam[] find_specific_families(FIGFam[] figfams, Strain[] strains) {
        ArrayList<FIGFam> fams = new ArrayList<>();
        for (FIGFam fam : figfams) {
            if (fam.hasStrain(strains, true)) {
                fams.add(fam);
            }
        }
        return fams.toArray(new FIGFam[]{});
    }

    private static FIGFam[] find_orthologues(FIGFam[] figfams, Strain[] strains, double ratioThreshold, boolean equals) {
        ArrayList<FIGFam> eqFams = new ArrayList<>();
        ArrayList<FIGFam> uneqFams = new ArrayList<>();
        int equal = 0, total = 0, exactly_one_perStrain = 0;
        for (FIGFam fam : figfams) {
            total++;
            if (fam.hasExactly_n_ProteinsPerMember(1)) {
                exactly_one_perStrain++;
                String p1 = fam.getMembers().get(strains[0]).toArray(new String[]{})[0];
                String p2 = fam.getMembers().get(strains[1]).toArray(new String[]{})[0];
                if (strains.length < 3) {
                    if (Diff.matchingRatio(proteome.get(p1), proteome.get(p2)) >= ratioThreshold) {
                        equal++;
                        eqFams.add(fam);
                    } else {
                        uneqFams.add(fam);
                    }
                } else {
                    String p3 = fam.getMembers().get(strains[2]).toArray(new String[]{})[0];
                    if ((Diff.matchingRatio(proteome.get(p1), proteome.get(p2)) >= ratioThreshold) && (Diff.matchingRatio(proteome.get(p1), proteome.get(p3)) >= ratioThreshold) && (Diff.matchingRatio(proteome.get(p2), proteome.get(p3)) >= ratioThreshold)) {
                        equal++;
                        eqFams.add(fam);
                    } else {
                        uneqFams.add(fam);
                    }
                }
            } else {
                uneqFams.add(fam);
            }
        }
        System.out.println("eq:" + equal + "\ntotal: " + total + " (" + exactly_one_perStrain + " with 1 protein per strain)");
        return equals ? eqFams.toArray(new FIGFam[]{}) : uneqFams.toArray(new FIGFam[]{});
    }

    public static void mafft(FIGFam[] fams) throws IOException, InterruptedException {
        int cnt = 0;
        for (int i = 0; i < fams.length; i++) {
            if (fams[i].getId().matches("FIG01306568|FIG00638284") || fams[i].hasExactly_n_ProteinsPerMember(1) || fams[i].getFunction().contains("ypothetical")) {continue;}
            PrintWriter pw = new PrintWriter("/tmp/test.fa");
            int length = 0;
            ArrayList<String> seq = new ArrayList<>();
            seqMap.put(fams [i].getId(), new ArrayList<String>());
            for (Map.Entry<Strain, HashSet<String>> entrySet : fams[i].getMembers().entrySet()) {
                HashSet<String> value = entrySet.getValue();
                for (String next : value) {
                    pw.write('>' + next + '\n');
                    seq.add(proteome.get(next));
                    seqMap.get(fams[i].getId()).add(next);
                    pw.write(proteome.get(next) + "\n\n");
                    length = Math.max(proteome.get(next).length(), length);//proteome.get(next).length() > length ? proteome.get(next).length() : length;
                }
            }
            pw.close();
            if(fams[i].getId().equalsIgnoreCase("FIG01366623")){
                    System.out.println("");
                }
            if (fams[i].hasExactly_n_ProteinsPerMember(1) && lcs_ratio(seq, length) < 0.8) {
                cnt = printRun_1_perFam(i, fams[i], length, cnt);
            }
            else if(!fams[i].hasExactly_n_ProteinsPerMember(1) && !lcs_moreThanOne_perFam(seq)){
                cnt = printRun(i, fams[i], length, cnt);
            }
        }
        System.out.println(cnt + " < 0.8");
    }
    
    private static ArrayList<Integer> clustersList = new ArrayList<>();// stores which sequence pairs are parallogous
    
    private static boolean lcs_moreThanOne_perFam(ArrayList<String> seq){
        ArrayList<HashSet<Integer>> clusters = new ArrayList<>();// stores which sequence pairs are parallogous
        clustersList = new ArrayList<>();
        HashSet<Integer> completeSet = new HashSet<>(), emptySet = new HashSet<>();
        for (int i = 0; i < seq.size(); i++) {
            completeSet.add(i);
        }
        for (int i = 0; i < seq.size(); i++) {
            for (int j = 0; j < seq.size(); j++) {
                if (j > i) {
                    String a = seq.get(i);
                    String b = seq.get(j);
                    String s = LCS.backtrack(LCS.fillMatrix(a, b), a, b, a.length(), b.length());
                    double lcs_ratio = 1.0 * s.length() / Math.max(a.length(), b.length());
                    if(lcs_ratio > 0.8){
                        boolean new_cluster = true;
                        for (HashSet<Integer> cluster : clusters) {
                            if(cluster.contains(i) || cluster.contains(j)){
                                new_cluster = false;
                                emptySet.add(i);
                                emptySet.add(j);
                                cluster.add(i);
                                cluster.add(j);
                                if(!clustersList.contains(i)){clustersList.add(i);}
                                if(!clustersList.contains(j)){clustersList.add(j);}
                            }
                        }
                        if(new_cluster){
                                HashSet<Integer> set = new HashSet<>();
                                set.add(i);
                                set.add(j);
                                emptySet.add(i);
                                emptySet.add(j);
                                clusters.add(set);
                                if(!clustersList.contains(i)){clustersList.add(i);}
                                if(!clustersList.contains(j)){clustersList.add(j);}
                            }
                    }
                }
            }
        }
        completeSet.removeAll(emptySet);
        for (Integer c : completeSet) {
            clustersList.add(c);
        }
        int members = -1;
        for (HashSet<Integer> cluster : clusters) {
            if(members == -1){
                members = cluster.size();
            }
            else if(members != cluster.size()){
                return false;
            }
        }
        return true;
    }
    
    private static double lcs_ratio(ArrayList<String> seq, int length){
        double lcs_ratio = 1.0;
            for (int j = 0; j < seq.size(); j++) {
                for (int k = 0; k < seq.size(); k++) {
                    if (k > j) {
                        String a = seq.get(j);
                        String b = seq.get(k);
                        lcs_ratio = Math.min(lcs_ratio, 1.0 * LCS.backtrack(LCS.fillMatrix(a, b), a, b, a.length(), b.length()).length() / length);
                    }
                }
            }
            return lcs_ratio;
    }

    private static int printRun(int i, FIGFam fam, int length, int cnt) throws IOException, InterruptedException {
        HashMap<String, StringBuilder> alignments = new HashMap<>();
        String[] sa = new String[]{"/home/tobias/mafft/bin/mafft", "--auto", "/tmp/test.fa"};
        System.out.println("Run " + i + " with " + fam.getId() + ": " + fam.getFunction());
        Process p = Runtime.getRuntime().exec(sa);
        p.waitFor();
        BufferedReader reader = new BufferedReader(new InputStreamReader(p.getInputStream()));
        StringBuilder sb = new StringBuilder().append('\n');
        String line;
        String id = "";
        while ((line = reader.readLine()) != null) {
            if (line.startsWith(">")) {
                id = line.substring(1, line.length());
                alignments.put(id, new StringBuilder());
            } else {
                alignments.get(id).append(line.trim());
            }
        }
        cnt++;
        for (Integer c : clustersList) {
            sb.append('>').append(seqMap.get(fam.getId()).get(c)).append('\n');
            sb.append(alignments.get(seqMap.get(fam.getId()).get(c))).append('\n');
        }
        sb.append("\n---\n");
        System.out.println(sb.toString());
        reader.close();
        return cnt;
    }
    
    private static int printRun_1_perFam(int i, FIGFam fam, int length, int cnt) throws IOException, InterruptedException {
        String[] sa = new String[]{"/home/tobias/mafft/bin/mafft", "--auto", "/tmp/test.fa"};
        System.out.println("Run " + i + " with " + fam.getId() + ": " + fam.getFunction());
        Process p = Runtime.getRuntime().exec(sa);
        p.waitFor();
        BufferedReader reader = new BufferedReader(new InputStreamReader(p.getInputStream()));
        StringBuilder sb = new StringBuilder().append('\n');
        String line;
        while ((line = reader.readLine()) != null) {
            if (line.startsWith(">")) {
                sb.append('\n').append(line).append('\n');
            } else {
                sb.append(line.trim());
            }
        }
        cnt++;
        System.out.println(sb.append("\n\n---\n\n").toString());
        reader.close();
        return cnt;
    }
    
    public static void main(String[] args) throws IOException, InterruptedException {
        String path = "/home/tobias/Dropbox/UNI/BACHELOR/Daten_Ergebnisse/";//"/home/h/harrert/Dropbox/UNI/BACHELOR/Daten_Ergebnisse/"
        read_proteome(path + "proteome/");
        FIGFam[] allFams = process_FigFam(path + "ProteinFamilyFeatures.txt");
        //Strain[] strains = new Strain[]{Strain.Q177, Strain.Q154};
        Strain[] strains = Strain.values();//test with all starins!
        FIGFam[] specificFams = find_specific_families(allFams, strains);
        FIGFam[] un_equalFams = find_orthologues(specificFams, strains, 0.98, false);//false means unequal
        mafft(un_equalFams);
        /*for (FIGFam equalFam : un_equalFams) {
         System.out.println(equalFam.getId());
         }*/
    }
}

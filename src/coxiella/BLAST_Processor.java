package coxiella;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;

/**
 *
 * @author harrert
 */
public class BLAST_Processor {

    private final HashMap<String, ProteinFunction> functionMap = new HashMap<>();
    private final double threshold = 1e-4;
    private final HashMap<String, String> equalProteins = new HashMap<>();
    private final HashMap<String, Boolean> equivProteinExists = new HashMap<>();
    
    
    /**
     * Check if two equal proteins exist in two different proteoms
     * <p>
     * Fills two maps with a bool map of all proteins and a map of all corresponding proteins.
     * @param file The input BLAST file
     * @throws IOException if file does not exist
     */
    public void proteinExists(String file) throws IOException{
        BufferedReader br = new BufferedReader(new FileReader(file));
        String line;
        boolean skip = false;
        int t=0,f=0, all=0, probe=0;
        while((line = br.readLine()) != null){
            boolean header = line.startsWith("#");
            String[] splt = line.split("\t");
            if(!header && splt[2].equals("100.00") && !splt[0].equals(splt[1])){
                String q = splt[0].split("@")[1];
                String s = splt[1].split("@")[1].substring(0, 1);
                if (!q.startsWith(s)) {
                    probe++;
                }
            }
            all = (line.startsWith("# Query:"))? all +1 : all;
            if(!header && !skip){
                String[] split = line.split("\t");
                String[] query = split[0].split("@");
                String[] subject = split[1].split("@");
                if(split[0].equals(split[1]) && !equivProteinExists.containsKey(query[1])){//self hit
                    equivProteinExists.put(query[1], false);
                    f++;
                }
                else if(!split[2].equals("100.00") && !equivProteinExists.containsKey(query[1])){ //no match
                    equivProteinExists.put(query[1], false);
                    skip = true;
                    f++;
                }
                else if(!query[0].equals(subject[0]) && split[2].equals("100.00")){//match found, no self hit
                    if(equivProteinExists.containsKey(query[1])){
                        f--;
                    }
                    equalProteins.put(query[1], subject[1]);
                    equivProteinExists.put(query[1], true);
                    skip = true;
                    t++;
                }
            }
            else if(header){
                skip = false;
            }
        }
        br.close();
        System.out.println("t: "+t+"\tf: "+f+"\tt+f: "+(t+f)+"\tall: "+all);
        System.out.println("Probe: "+probe);
    }
    
    /**
     * Process a BLASTp file
     * <p>
     * Assign function(s) to proteinIDs from the BLAST results
     * @param file The path to the BLAST file
     * @throws IOException if file does not exist etc.
     */
    public void processBLAST(String file) throws IOException{
        BufferedReader br = new BufferedReader(new FileReader(file));
        String line, proteinID = "", function="";
        while((line = br.readLine()) != null){
            if(line.startsWith("# Query:")){
                function = line.split("@\\d+\\s+")[1];
                if(!function.contains("hypothetical protein")){
                    proteinID = line.split(" ")[2].split("@")[1];
                    functionMap.put(proteinID, new ProteinFunction(proteinID, function));
                }
            }
            else if(!line.startsWith("#") && !function.contains("hypothetical protein")){
                final String[] tabs = line.split("\t");
                final String target = tabs[1].split("@")[1];
                if(!proteinID.equals(target) && Double.compare(Double.parseDouble(tabs[10]), threshold) <= 0){
                    if(functionMap.containsKey(target) && !functionMap.get(target).getFunction().contains(function)){
                        functionMap.get(target).addFunction(function);
                    }
                    else {
                        functionMap.put(target, new ProteinFunction(proteinID, function));
                    }
                }
            }
        }
        br.close();
    }
    
    /**
     * Prints the contents of the map containing the proteinIDs' function(s) to a file
     * <p>
     * The format is tab separated: ID \t number of functions (#) \t function(s)
     * @param file The path to the output file or where it will be created
     * @throws FileNotFoundException 
     */
    public void printProteinFunction(String file) throws FileNotFoundException{
        StringBuilder sb = new StringBuilder("ID\t#\tfunction\n");
        for (Map.Entry<String, ProteinFunction> entrySet : functionMap.entrySet()) {
            String id = entrySet.getKey();
            ProteinFunction value = entrySet.getValue();
            sb.append(id);
            ArrayList<String> function = value.getFunction();
            sb.append('\t').append(function.size());
            for (String s : function) {
                sb.append('\t').append(s);
            }
            sb.append('\n');
        }
        PrintWriter pw = new PrintWriter(file);
        pw.write(sb.toString());
        pw.close();
    }
    
    public static void main(String[] args) throws IOException {
        BLAST_Processor bp = new BLAST_Processor();
        bp.proteinExists("/Users/Tobias/ThreeByThree.csv");///home/h/harrert/Desktop/Coxiella/BLAST_out/ThreeByThree.tab
        //bp.processBLAST("/home/h/harrert/Desktop/Coxiella/BLAST_out/ThreeByThree.tab");
        //bp.printProteinFunction("/home/h/harrert/Desktop/cox.csv");
    }
}

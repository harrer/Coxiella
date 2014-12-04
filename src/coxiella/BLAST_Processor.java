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
        while((line = br.readLine()) != null){
            boolean header = line.startsWith("#");
            if(!header && !skip){
                String[] split = line.split("\t");
                String[] query = split[0].split("@");
                String[] subject = split[1].split("@");
                if(!split[2].equals("100.00")){
                    equivProteinExists.put(query[1], false);
                    skip = true;
                }
                else if(!query[0].equals(subject[0])){//match found, no self hit
                    equalProteins.put(query[1], subject[1]);
                    equivProteinExists.put(query[1], true);
                    skip = true;
                }
            }
            else if(header){
                skip = false;
            }
        }
        br.close();
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
        bp.proteinExists("/Users/Tobias/ThreeByThree.csv");
        //bp.processBLAST("/home/h/harrert/Desktop/Coxiella/BLAST_out/ThreeByThree.tab");
        //bp.printProteinFunction("/home/h/harrert/Desktop/cox.csv");
    }
}

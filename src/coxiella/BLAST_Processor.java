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
        bp.processBLAST("/home/h/harrert/Desktop/Coxiella/BLAST_out/ThreeByThree.tab");
        bp.printProteinFunction("/home/h/harrert/Desktop/cox.csv");
    }
}

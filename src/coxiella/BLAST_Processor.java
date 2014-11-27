package coxiella;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashMap;

/**
 *
 * @author harrert
 */
public class BLAST_Processor {

    private final HashMap<String, ProteinFunction> functionMap = new HashMap<>();
    private final double threshold = 1e-4;
    
    public void processBLAST(String file) throws IOException{
        BufferedReader br = new BufferedReader(new FileReader(file));
        String line, proteinID = "", function="";
        while((line = br.readLine()) != null){
            if(line.startsWith("# Query:")){
                function = line.split("@\\d+\\s+")[1];
                proteinID = line.split(" ")[2].split("@")[1];
                functionMap.put(proteinID, new ProteinFunction(proteinID, function));
            }
            else if(!line.startsWith("#") && !function.equalsIgnoreCase("hyp")){//////////////////
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
    
    public static void main(String[] args) throws IOException {
        BLAST_Processor bp = new BLAST_Processor();
        bp.processBLAST("/home/h/harrert/Desktop/Coxiella/BLAST_out/ThreeByThree.tab");
        System.out.println("");
    }
}

package coxiella;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.HashMap;
import java.util.Map;

/**
 *
 * @author Tobias
 */
public class Coxiella {
    
    private final HashMap<String, String> seq_map = new HashMap<>(), fun_map = new HashMap<>();
    private final HashMap<String, Integer> functionCountMap = new HashMap<>();
    
    public void trimFastaHeader(String file) throws IOException{
        BufferedReader br = new BufferedReader(new FileReader(file));
        StringBuilder sb = new StringBuilder();
        String line;
        while((line = br.readLine()) != null){
            if(line.startsWith(">")){//fasta header
                sb.append(line.split("@")[1]);
            }
            else{//sequence
                sb.append(line);
            }
        }
        PrintWriter pw = new PrintWriter(file);
        pw.write(sb.toString());
    }
    
    public void extractFasta(String file) throws IOException{
        BufferedReader br = new BufferedReader(new FileReader(file));
        StringBuilder sb = new StringBuilder();
        String line, key, function;
        
        Integer i = 0;///////////////////
        line = br.readLine();
        key = line.split("@")[1];
        function = line.split("   ")[1];
        fun_map.put(key, function);
        while((line = br.readLine()) != null){
            if(line.startsWith(">")){
                seq_map.put(key, sb.toString());
                sb.delete(0, sb.length());
                key = line.split("@")[1];
                function = line.split("   ")[1];
                fun_map.put(key, function);
                i = functionCountMap.containsKey(function)? functionCountMap.get(function) + 1 : 1;//////////////
                functionCountMap.put(function, i);//////////////////
            }
            else{
                sb.append(line);
            }
        }
        seq_map.put(key, sb.toString());
        br.close();
        i=0; int j=0;
        HashMap<Integer, Integer> m = new HashMap<>();
        for (Map.Entry<String, Integer> entrySet : functionCountMap.entrySet()) {
            key = entrySet.getKey();
            Integer value = entrySet.getValue();
            if(value > 1){
                i++;
                System.out.println(key + '\t' + value + '\n');
                j = m.containsKey(value)? m.get(value) + 1 : 1;
                m.put(value, j);
            }
        }
        System.out.println(i + " entries");
        System.out.println("Distribution:\n");
        for (Map.Entry<Integer, Integer> entrySet : m.entrySet()) {
            Integer key1 = entrySet.getKey();
            Integer value = entrySet.getValue();
            System.out.println(key1 + "\t" + value);
        }
    }

    /**
     * @param args the command line arguments
     * @throws java.io.IOException
     */
    public static void main(String[] args) throws IOException {
        // TODO code application logic here
        Coxiella cox = new Coxiella();
        cox.extractFasta("/home/h/harrert/Desktop/p3_r57631_Cox_burne.fa");//"/Users/Tobias/Desktop/p3_r57631_Cox_burne.fa"
    }
    
}

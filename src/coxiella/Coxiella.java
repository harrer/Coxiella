package coxiella;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashMap;

/**
 *
 * @author Tobias
 */
public class Coxiella {
    
    public void extractFasta(String file) throws IOException{
        BufferedReader br = new BufferedReader(new FileReader(file));
        StringBuilder sb = new StringBuilder();
        String line, key, function;
        HashMap<String, String> seq_map = new HashMap<>(), fun_map = new HashMap<>();
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
            }
            else{
                sb.append(line);
            }
        }
        seq_map.put(key, sb.toString());
        br.close();
    }

    /**
     * @param args the command line arguments
     * @throws java.io.IOException
     */
    public static void main(String[] args) throws IOException {
        // TODO code application logic here
        Coxiella cox = new Coxiella();
        cox.extractFasta("/home/tobias/Dropbox/UNI/BACHELOR/Sequenzen/RSA_493_NM1/p3_r57631_Cox_burne.fa");//("/Users/Tobias/Desktop/p3_r57631_Cox_burne.fa");
    }
    
}

package coxiella;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
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
        String line, key, sequence;
        HashMap<String, String> map = new HashMap<>();
        key = br.readLine().split("@")[1];
        while((line = br.readLine()) != null){
            if(line.startsWith(">")){
                map.put(key, sb.toString());
                sb.delete(0, sb.length());
                key = line.split("@")[1];
            }
            else{
                sb.append(line);
            }
        }
        br.close();
    }

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) throws IOException {
        // TODO code application logic here
        Coxiella cox = new Coxiella();
        cox.extractFasta("/Users/Tobias/Downloads/p3_r57631_Cox_burne.fa");
        System.out.println("");
    }
    
}

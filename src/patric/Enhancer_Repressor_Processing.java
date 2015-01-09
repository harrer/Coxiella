package patric;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;

/**
 *
 * @author tobias
 */
public class Enhancer_Repressor_Processing {
    
    private int[] lineToInt(String line){
        String [] prot = line.split("\t")[5].split(" ")[1].split("\\(")[1].split("\\)")[0].split("\\|");
        return new int[]{Integer.parseInt(prot[0]), Integer.parseInt(prot[1]), Integer.parseInt(prot[2])};
    }
    
    public void processFile(String file, String outPath, boolean ommit_noFun) throws IOException{
        BufferedReader br = new BufferedReader(new FileReader(file));
        br.readLine();
        StringBuilder enhancer = new StringBuilder("FigFam\tQ177|Q154|RSA_493\tPatric Locus tags Q177 (107188_...)\tQ154 (77120_...)\tRSA_493 (82552_...)\tnumber of proteins\tortholog\tfunction\n");
        StringBuilder repressor = new StringBuilder("FigFam\tQ177|Q154|RSA_493\tPatric Locus tags Q177 (107188_...)\tQ154 (77120_...)\tRSA_493 (82552_...)\tnumber of proteins\tortholog\tfunction\n");
        String line;
        while((line = br.readLine()) != null && !line.equals("")){
            int[] n_of_proteins = lineToInt(line); //stores content of "number of proteins": eg 10(4|4|2) -> 4,4,2
            String[] columns = line.split("\t");
            if(n_of_proteins[0] > n_of_proteins[2] || n_of_proteins[1] > n_of_proteins[2] || columns[1].equalsIgnoreCase("100") || columns[1].equalsIgnoreCase("010")){
                if(ommit_noFun){
                    if(!columns[7].equalsIgnoreCase("---")){
                        repressor.append(line).append('\n');
                    }
                }
                else{
                    repressor.append(line).append('\n');
                }
            }
            else if((n_of_proteins[0] < n_of_proteins[2] && n_of_proteins[1] < n_of_proteins[2]) || columns[1].equalsIgnoreCase("001")){
                if(ommit_noFun){
                    if(!columns[7].equalsIgnoreCase("---")){
                        enhancer.append(line).append('\n');
                    }
                }
                else{
                    enhancer.append(line).append('\n');
                }
            }
            else{
                System.out.println(line);
            }
        }
        br.close();
        PrintWriter pw = new PrintWriter(outPath+"enhancers.csv");
        pw.write(enhancer.toString());
        pw.close();
        pw = new PrintWriter(outPath+"repressors.csv");
        pw.write(repressor.toString());
        pw.close();
    }
    
    public static void main(String[] args) throws IOException {
        String dir = "/home/tobias/Desktop/";
        new Enhancer_Repressor_Processing().processFile(dir+"figfam_genomes.csv", dir, true);
    }
}

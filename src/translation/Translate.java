/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package translation;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashMap;

/**
 *
 * @author harrert
 */
public class Translate {
    
    private final HashMap <String, String> h=new HashMap <>();
    
    public Translate(boolean shortCode){
        if(shortCode){init_h_short();}
        else{init_h();}
    }
    
    private void init_h_short(){
        h.put("UUU", "F");
        h.put("UUC", "F");
        h.put("UUA", "L");
        h.put("UUG", "<L>");
        h.put("CUU", "L");
        h.put("CUC", "L");
        h.put("CUA", "L");
        h.put("CUG", "<L>");
        h.put("AUU", "<I>");
        h.put("AUC", "<I>");
        h.put("AUA", "<I>");
        h.put("AUG", "<M>");
        h.put("GUU", "V");
        h.put("GUC", "V");
        h.put("GUA", "V");
        h.put("GUG", "<V>");
        h.put("UCU", "S");
        h.put("UCC", "S");
        h.put("UCA", "S");
        h.put("UCG", "S");
        h.put("CCU", "P");
        h.put("CCC", "P");
        h.put("CCA", "P");
        h.put("CCG", "P");
        h.put("ACU", "T");
        h.put("ACC", "T");
        h.put("ACA", "T");
        h.put("ACG", "T");
        h.put("GCU", "A");
        h.put("GCC", "A");
        h.put("GCA", "A");
        h.put("GCG", "A");
        h.put("UAU", "Y");
        h.put("UAC", "Y");
        h.put("UAA", "<Stop>");
        h.put("UAG", "<Stop>");
        h.put("CAU", "H");
        h.put("CAC", "H");
        h.put("CAA", "Q");
        h.put("CAG", "Q");
        h.put("AAU", "N");
        h.put("AAC", "N");
        h.put("AAA", "K");
        h.put("AAG", "K");
        h.put("GAU", "D");
        h.put("GAC", "D");
        h.put("GAA", "E");
        h.put("GAG", "E");
        h.put("UGU", "C");
        h.put("UGC", "C");
        h.put("UGA", "<Stop>");
        h.put("UGG", "W");
        h.put("CGU", "R");
        h.put("CGG", "R");//hat gefehlt bzw. war CGC; wurde am 09.11.10 als CGG richtiggestellt, nachdem bei der Translation eine NullPointerException aufgetaucht war.
        h.put("CGA", "R");
        h.put("CGC", "R");
        h.put("AGU", "S");
        h.put("AGC", "S");
        h.put("AGA", "R");
        h.put("AGG", "R");
        h.put("GGU", "G");
        h.put("GGC", "G");
        h.put("GGA", "G");
        h.put("GGG", "G");
    }
    
    private void init_h(){
        h.put("UUU", "Phe");
        h.put("UUC", "Phe");
        h.put("UUA", "Leu");
        h.put("UUG", "<Leu>");
        h.put("CUU", "Leu");
        h.put("CUC", "Leu");
        h.put("CUA", "Leu");
        h.put("CUG", "<Leu>");
        h.put("AUU", "<Ile>");
        h.put("AUC", "<Ile>");
        h.put("AUA", "<Ile>");
        h.put("AUG", "<Met>");
        h.put("GUU", "Val");
        h.put("GUC", "Val");
        h.put("GUA", "Val");
        h.put("GUG", "<Val>");
        h.put("UCU", "Ser");
        h.put("UCC", "Ser");
        h.put("UCA", "Ser");
        h.put("UCG", "Ser");
        h.put("CCU", "Pro");
        h.put("CCC", "Pro");
        h.put("CCA", "Pro");
        h.put("CCG", "Pro");
        h.put("ACU", "Thr");
        h.put("ACC", "Thr");
        h.put("ACA", "Thr");
        h.put("ACG", "Thr");
        h.put("GCU", "Ala");
        h.put("GCC", "Ala");
        h.put("GCA", "Ala");
        h.put("GCG", "Ala");
        h.put("UAU", "Tyr");
        h.put("UAC", "Tyr");
        h.put("UAA", "<Stop>");
        h.put("UAG", "<Stop>");
        h.put("CAU", "His");
        h.put("CAC", "His");
        h.put("CAA", "Gln");
        h.put("CAG", "Gln");
        h.put("AAU", "Asn");
        h.put("AAC", "Asn");
        h.put("AAA", "Lys");
        h.put("AAG", "Lys");
        h.put("GAU", "Asp");
        h.put("GAC", "Asp");
        h.put("GAA", "Glu");
        h.put("GAG", "Glu");
        h.put("UGU", "Cys");
        h.put("UGC", "Cys");
        h.put("UGA", "<Stop>");
        h.put("UGG", "Trp");
        h.put("CGU", "Arg");
        h.put("CGG", "Arg");//hat gefehlt bzw. war CGC; wurde am 09.11.10 als CGG richtiggestellt, nachdem bei der Translation eine NullPointerException aufgetaucht war.
        h.put("CGA", "Arg");
        h.put("CGC", "Arg");
        h.put("AGU", "Ser");
        h.put("AGC", "Ser");
        h.put("AGA", "Arg");
        h.put("AGG", "Arg");
        h.put("GGU", "Gly");
        h.put("GGC", "Gly");
        h.put("GGA", "Gly");
        h.put("GGG", "Gly");
    }
    
    public String translation (String mRNA) {
        mRNA = mRNA.replaceAll("\n", "");
        mRNA = to_mRNA(mRNA);
        StringBuilder pr = new StringBuilder();
        for (int i = 0; i < 3; i++) {
            pr.append("frame ").append(i+1).append(":\n");
            for (int j = i; j < mRNA.length(); j += 3) {
                if (j < mRNA.length() - 2) {
                    pr.append(h.get(mRNA.substring(j, j + 3)));//.append(' ');
                } else {
                    break;
                }
            }
            pr.append("\n\n");
        }
        //complementary strand:
        mRNA = compStrand(mRNA);
        for (int i = 0; i < 3; i++) {
            pr.append("complementary frame ").append(i+4).append(":\n");
            for (int j = i; j < mRNA.length(); j += 3) {
                if (j < mRNA.length() - 2) {
                    pr.append(h.get(mRNA.substring(j, j + 3)));//.append(' ');
                } else {
                    break;
                }
            }
            pr.append("\n\n");
        }
        return pr.toString();
    }
    
    private String compStrand(String strand){
        strand = new StringBuilder().append(strand).reverse().toString();
        StringBuilder sb = new StringBuilder();
        for (int i = 0; i < strand.length(); i++) {
            switch (strand.charAt(i)) {
                case 'A':
                    sb.append('U');
                    break;
                case 'U':
                    sb.append('A');
                    break;
                case 'G':
                    sb.append('C');
                    break;
                case 'C':
                    sb.append('G');
                    break;
            }
        }
        return sb.toString();
    }
    
    private String to_mRNA(String dna){
        dna = dna.toUpperCase();
        StringBuilder sb = new StringBuilder();
        for (int i = 0; i < dna.length(); i++) {
            sb.append(dna.charAt(i) == 'T'? 'U' : dna.charAt(i));
        }
        return sb.toString();
    }
    
    public static String readFromFile(String file) throws IOException{
        BufferedReader br = new BufferedReader(new FileReader(file));
        String line;
        StringBuilder sb = new StringBuilder();
        while((line = br.readLine()) != null){
            sb.append(line);
        }
        return sb.toString();
    }
    
    public static void main(String[] args) throws IOException {
        //String s = "atggaaatcttctattttgttgccacattgaaaagtttttcccgggccgctttagaactgggtgtttcgaaaggatatgtcagcactcaaattaccgcattagaaaaaga";
        String s = "AATGTTAGGGCATTAG";
        String c = "TTACAATCCCGTAATC";
        Translate t = new Translate(true);
        System.out.println(t.translation(readFromFile("/tmp/test.txt")));
    }
}

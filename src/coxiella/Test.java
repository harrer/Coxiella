/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package coxiella;

/**
 *
 * @author harrert
 */
public class Test {
    
    public static void main(String[] args) {
        String line = "FIG00551682     110     107188_1198     77120_0362      ---     2 (1|1|0)       1       ---";
        String [] prot = line.split("\t")[5].split(" ")[1].split("\\(")[1].split("\\)")[0].split("|");
        //return new int[]{Integer.parseInt(prot[0]), Integer.parseInt(prot[1]), Integer.parseInt(prot[2])};
        //System.out.println("9 (4|4|1)".split(" ")[1].split("\\(")[1].split("\\)")[0]);
    }
}

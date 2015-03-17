package lcs;

/**
 *
 * @author tobias
 */
public class LCS {
    
    public static int[][] fillMatrix(String a, String b){
        int[][] matrix = new int[a.length()+1][b.length()+1];
        for (int i = 1; i <= a.length(); i++) {
            for (int j = 1; j <= b.length(); j++) {
                matrix[i][j] = (a.charAt(i-1) == b.charAt(j-1))? matrix[i-1][j-1] + 1 : Math.max(matrix[i][j-1], matrix[i-1][j]);
            }
        }
        return matrix;
    }
    
    public static String backtrack(int[][] matrix, String a, String b, int i, int j){
        return (i == 0 || j == 0)? "" : (a.charAt(i-1) == b.charAt(j-1))? backtrack(matrix, a, b, i-1, j-1) + a.charAt(i-1) : (matrix[i][j-1] > matrix[i-1][j])? backtrack(matrix, a, b, i, j-1) : backtrack(matrix, a, b, i-1, j);
    }
    
    public static void main(String[] args) {
        String a = "XMJYAUZ";
        String b = "MZJAWXU";
        int[][] m = fillMatrix(a, b);
        String s = backtrack(m, a, b, a.length(), b.length());
        System.out.println(1.0*s.length()/Math.max(a.length(), b.length()));
    }
    
}

package translation;

/**
 *
 * @author tobias
 */
public class Diff {
    
    public static String showDiffs(String a, String b, boolean invert){
        StringBuilder sb = new StringBuilder(a).append('\n');
        char c = invert? ' ' : '|';
        char nc = invert? '|' : ' ';
        for (int i = 0; i < Math.min(a.length(),b.length()); i++) {
            sb.append(a.charAt(i) == b.charAt(i)? c : nc);
        }
        return sb.append('\n').append(b).toString();
    }
    
    public static double matchingRatio(String a, String b){
        double matching = 0;
        int length = Math.min(a.length(),b.length());
        for (int i = 0; i < length; i++) {
            if(a.charAt(i) == b.charAt(i)){
                matching++;
            }
        }
        return matching/length;
    }
    
    public static void main(String[] args) {
        String a = "MEAQSSHQAPLSPKGSLGLFRGYISILRELFTILSEWCSPTTMWMGLVLPGIFINLLSLAFPLILLQVYDRVIPQQSIYTLTFLIVGGIIVTGIAMVLSILRSVSINWTSARFEYFTHLRVFQRLLACKLDDFQKEGSGYY";
        String b = "MEAQSSHQAPLSPKGSLGLFRGYINILRELFTILSEWCSPTTMWMGLVLPGIFINLLSLAFPLILLQVYDRVIPQQSIYTLTFLIVGGIIVTGIAMVLSILRSVSINWTSARFEYFTHLRVFQRLLACKLDDFQKEGSG";
        System.out.println(showDiffs(a, b, true));
    }
}

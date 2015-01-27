package patric;

import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;

/**
 *
 * @author tobias
 */
public class FIGFam {
    
    private final String id;
    private final HashMap<Strain, HashSet<String>> members = new HashMap<>();

    public FIGFam(String id) {
        this.id = id;
    }

    public String getId() {
        return id;
    }

    public HashMap<Strain, HashSet<String>> getMembers() {
        return members;
    }
    
    public void addMember(String protein, Strain strain){
        HashSet proteins = members.getOrDefault(strain, new HashSet<String>());
        proteins.add(protein);
        members.put(strain, proteins);
    }
    
    public boolean hasStrain(Strain strain){
        return members.containsKey(strain);
    }
    
    public boolean hasStrain(Strain[] strains, boolean exclusively){
        if (exclusively) {
            return members.keySet().equals(new HashSet(Arrays.asList(strains)));
        } else {
            for (Strain strain : strains) {
                if (!hasStrain(strain)) {
                    return false;
                }
            }
            return true;
        }
    }
    
    public boolean hasExactly_n_ProteinsPerMember(int n){
        for (Map.Entry<Strain, HashSet<String>> entrySet : members.entrySet()) {
            HashSet<String> value = entrySet.getValue();
            if(value.size() != n){
                return false;
            }
        }
        return true;
    }
    
}

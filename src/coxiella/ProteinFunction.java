package coxiella;

import java.util.ArrayList;

/**
 *
 * @author harrert
 */
public class ProteinFunction {
    
    private final String id;
    private final ArrayList<String> function = new ArrayList<>();
    
    public ProteinFunction(String proteinID, String function){
        id = proteinID;
        this.function.add(function);
    }
    
    public boolean equals(ProteinFunction proteinFunction){
        return proteinFunction.id.equals(this.id);
    }

    public String getId() {
        return id;
    }

    public boolean addFunction(String function){
        if(this.function.contains(function)){
            return false;
        }
        else{
            this.function.add(function);
            return true;
        }
    }

    public ArrayList<String> getFunction() {
        return function;
    }
    
    
}

/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package refineassemblydraft;

import datastruct.AlignmentLine;
import java.util.Comparator;

/**
 *
 * @author desktop
 */
public class SortScaffoldByScaffoldCoords implements Comparator<AlignmentLine>{

    @Override
    public int compare(AlignmentLine t, AlignmentLine t1) {
        if(t1.Name().equals(t.Name())){
            if(t.reversedCoords() && t1.reversedCoords()){
                return t1.ContigEnd() - t.ContigEnd();
            }else if(t.reversedCoords() && !t1.reversedCoords()){
                return t1.ContigStart() - t.ContigEnd();
            }else if(!t.reversedCoords() && t1.reversedCoords()){
                return t1.ContigEnd() - t.ContigStart();
            }else{
                return t1.ContigStart() - t.ContigStart();
            }
        }else{
            return t1.Name().compareTo(t.Name());
        }
    }
    
}

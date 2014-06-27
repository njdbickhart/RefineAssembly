/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package datastruct;

import factories.AlignmentLineFactory;
import file.BedAbstract;
import file.BedMap;
import java.util.ArrayList;

/**
 *
 * @author bickhart
 */
public class AGPChrMap extends BedMap{
    
    public ArrayList<AlignmentLine> getSortedAlignments(String chr){
        ArrayList<AlignmentLine> lines = new ArrayList<>();
        for(BedAbstract b : this.getSortedBedAbstractList(chr)){
            lines.add((AlignmentLine) b);
        }
        return lines;
    }
    
    public void addAlignmentLines(AlignmentLineFactory fac){
        for(String contig : fac.getContigList()){
            for(AlignmentLine l : fac.collectAlignFromContig(contig)){
                this.addBedData(l);
            }
        }
    }
}

/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package datastruct;

import file.BedMap;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import refineassemblydraft.SortScaffoldByScaffoldCoords;

/**
 *
 * @author desktop
 */
public class SuperScaffold implements Comparable<SuperScaffold>{
    private ArrayList<String> scaffNames = new ArrayList<>();
    private ArrayList<AlignmentLine> alignments = new ArrayList<>();
    public int length = 0;
    
    public void addAlignments(ArrayList<AlignmentLine> aligns, ArrayList<String> names){
        this.scaffNames.addAll(names);
        this.alignments.addAll(aligns);
    }
    
    public void addAlignments(AlignmentLine align, String name){
        
            this.scaffNames.add(name);
            this.alignments.add(align);
    }
   
    
    public ArrayList<AlignmentLine> getOutScaffoldListSegments(){
        //Collections.sort(this.alignments, new SortScaffoldByScaffoldCoords());
        return this.alignments;
    }
    
    public boolean overlaps(AlignmentLine l, int ovlpBaseNum, HashSet<String> genuines){
        for(int x = 0; x < alignments.size(); x++){
            AlignmentLine t = alignments.get(x);
            boolean missStart = ((l.missStart && t.missEnd) || (scaffoldUMDBeginning(l) || scaffoldUMDEnd(t)))? true : false;
            boolean missEnd = ((l.missEnd && t.missStart) || (scaffoldUMDEnd(l) || scaffoldUMDBeginning(t)))? true : false;
            
            if(l.Name().equals("7180008669549")){
                //System.err.println("hey");
            }
            
            int ovlp = overlap(t.Start(), t.End(), l.Start(), l.End());
            if(t.Chr().equals(l.Chr())){
                if(!t.Name().equals(l.Name()) &&
                       ((missStart && l.Start() > t.Start()) || (missEnd && l.Start() < t.End() )) &&
                        ovlp >= -ovlpBaseNum){
                    // add logic to add this alignment line to the current list in appropriate order and 
                    // change the SuperScaffoldBuilder to avoid duplicating the addition!
                    if(l.Name().equals("7180008668866")){
                        //System.err.println("hey");
                    }
                    if(missEnd && l.Start() < t.Start()){
                        alignments.add(x, l);
                        scaffNames.add(x, l.Name());
                    }else{
                        alignments.add(x + 1, l);
                        scaffNames.add(x + 1, l.Name());
                    }
                        
                    return true;
                }else if(t.Name().equals(l.Name()) &&
                        ((l.Start() > t.Start() && !l.missStart && !t.missEnd) ||
                        (l.Start() < t.Start() && !t.missEnd && !t.missStart)) &&
                        (overlap(l.contigStart, l.contigEnd, t.contigStart, t.contigEnd, l.reversedCoords(), t.reversedCoords()) >= -ovlpBaseNum ||
                        (genuineMatchStartA(l,t) || genuineMatchEndA(l,t)))){
                    
                    if(genuineMatchEndA(l,t) && earlierScaffold(l,t)){
                        alignments.add(x, l);
                        scaffNames.add(x, l.Name());
                    }else if(genuineMatchStartA(l,t) && !earlierScaffold(l,t)){
                        alignments.add(x + 1, l);
                        scaffNames.add(x + 1, l.Name());
                    }else{
                        return false;
                    }
                    return true;
                }
            }else if(t.Name().equals(l.Name())){
                if(l.Name().equals("7180008669549")){
                    System.err.println("hey");
                }
                if((genuineMatchStartA(l,t) || genuineMatchEndA(l,t)) &&
                    overlap(l.contigStart, l.contigEnd, t.contigStart, t.contigEnd, l.reversedCoords(), t.reversedCoords()) >= -ovlpBaseNum){
                    
                    //Check if this alignment is against a larger alignment list
                    if(this.alignments.size() > 1){                        
                        if(x == 0 && x != this.alignments.size() - 1){
                            AlignmentLine next = alignments.get(x + 1);
                            if(t.Chr().equals(next.Chr()) && t.Start() < next.Start()){
                                alignments.add(x, l);
                                scaffNames.add(x, l.Name());
                            }else{
                                alignments.add(x + 1, l);
                                scaffNames.add(x + 1, l.Name());
                            }
                        }else if(x > 0 && x != this.alignments.size() - 1){
                            AlignmentLine next = alignments.get(x + 1);
                            AlignmentLine previous = alignments.get(x - 1);
                            if(t.Chr().equals(next.Chr()) && t.Start() < next.Start() && 
                                    !t.Chr().equals(previous.Chr())){
                                alignments.add(x, l);
                                scaffNames.add(x, l.Name());
                            }else if(t.Chr().equals(previous.Chr()) && t.Start() < previous.Start() && 
                                    !t.Chr().equals(next.Chr())){
                                alignments.add(x + 1, l);
                                scaffNames.add(x + 1, l.Name());
                            }else{
                                // We cannot accurately place this one, so let's put it all the way at the end
                                alignments.add(l);
                                scaffNames.add(l.Name());
                            }
                        }else{
                            AlignmentLine previous = alignments.get(x - 1);
                            if(t.Chr().equals(previous.Chr()) && t.Start() < previous.Start()){
                                alignments.add(x, l);
                                scaffNames.add(x, l.Name());
                            }else{
                                alignments.add(x + 1, l);
                                scaffNames.add(x + 1, l.Name());
                            }
                        }
                    }else{
                        // Just a small list; going to sort by scaffold alignment order
                        if(genuineMatchEndA(l,t) || earlierScaffold(l, t)){
                            alignments.add(x, l);
                            scaffNames.add(x, l.Name());
                        }else{
                            alignments.add(x + 1, l);
                            scaffNames.add(x + 1, l.Name());
                        }
                    }
                    return true;
                }
            }
        }
        return false;
    }
    
    private boolean scaffoldUMDBeginning(AlignmentLine a){
        if((a.reversedCoords() && a.contigSize - a.contigEnd < 5000)
                || (!a.reversedCoords() && a.contigStart < 5000))
            return true;
        return false;
    }
    
    private boolean scaffoldUMDEnd(AlignmentLine a){
        if((a.reversedCoords() && a.contigStart < 5000)
                || (!a.reversedCoords() && a.contigSize - a.contigEnd < 5000))
            return true;
        return false;
    }
    
    private boolean earlierScaffold(AlignmentLine a, AlignmentLine b){
        if(a.reversedCoords() && b.reversedCoords() && a.contigEnd < b.contigEnd){
            return true;
        }else if(a.reversedCoords() && !b.reversedCoords() && a.contigEnd < b.contigStart){
            return true;
        }else if(!a.reversedCoords() && b.reversedCoords() && a.contigStart < b.contigEnd){
            return true;
        }else if(!a.reversedCoords() && !b.reversedCoords() && a.contigStart < b.contigStart){
            return true;
        }
        return false;
    }
    
    private boolean genuineMatchStartA(AlignmentLine a, AlignmentLine b){
        if(a.reversedCoords() && b.reversedCoords() && a.contigEnd > b.contigEnd && a.genEnd && b.genStart){
            return true;
        }else if(a.reversedCoords() && !b.reversedCoords() && a.contigEnd > b.contigStart && a.genEnd && b.genEnd){
            return true;
        }else if(!a.reversedCoords() && b.reversedCoords() && a.contigStart > b.contigEnd && a.genStart && b.genStart){
            return true;
        }else if(!a.reversedCoords() && !b.reversedCoords() && a.contigStart > b.contigStart && a.genStart && b.genEnd){
            return true;
        }
        return false;
    }
    
    private boolean genuineMatchEndA(AlignmentLine a, AlignmentLine b){
        if(a.reversedCoords() && b.reversedCoords() && a.contigEnd < b.contigEnd && a.genStart && b.genEnd){
            return true;
        }else if(a.reversedCoords() && !b.reversedCoords() && a.contigEnd < b.contigStart && a.genStart && b.genStart){
            return true;
        }else if(!a.reversedCoords() && b.reversedCoords() && a.contigStart < b.contigEnd && a.genEnd && b.genEnd){
            return true;
        }else if(!a.reversedCoords() && !b.reversedCoords() && a.contigStart < b.contigStart && a.genEnd && b.genStart){
            return true;
        }
        return false;
    }
    
    public boolean containsAlign(AlignmentLine a){
        for(AlignmentLine l : this.alignments){
            if(l.equals(a))
                return true;
        }
        return false;
    }
    
    public void calculateLength(){
        int len = 0;
        HashMap<String, Integer> names = new HashMap<>();
        for(AlignmentLine l : alignments){
            names.put(l.Name(), l.contigSize);
        }
        
        for(String name : names.keySet()){
            len += names.get(name);
        }
        this.length = len;
    }
    
    @Override
    public int compareTo(SuperScaffold t) {
        if(t.length == 0){
            t.calculateLength();
        }
        if(this.length == 0){
            this.calculateLength();
        }
        return t.length - this.length;
    }
    
    private int least(int a, int b){
        if(a < b) {
            return a;
        }
        else {
            return b;
        }
    }
    private int most(int a, int b){
        if(a > b) {
            return a;
        }
        else {
            return b;
        }
    }
    private int overlap(int s1, int e1, int s2, int e2, boolean rev1, boolean rev2){
        int ts1 = s1;
        int te1 = e1;
        int ts2 = s2;
        int te2 = e2;
        if(rev1){
            ts1 = e1;
            te1 = s1;
        }
        if(rev2){
            ts2 = e2;
            te2 = s2;
        }
        return (least(te1, te2) - most(ts1, ts2));
    }
    private int overlap(int s1, int e1, int s2, int e2){
        return (least(e1, e2) - most(s1, s2));
    }
}

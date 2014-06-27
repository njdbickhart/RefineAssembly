/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package datastruct;

import file.BedAbstract;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 *
 * @author Dbick
 */
public class AlignmentLine extends BedAbstract{
    protected int contigStart;
    protected int contigEnd;
    protected int contigSize;
    //Using the "Name" field for the scaffold Name
    protected float percGenomeAlign; //Proportion of the scaffold that aligns to the genome
    protected float percSegmentAlign; //Alignment of the relevant scaffold segment to the genome
    protected int proportionSegAlign; //The BP size difference between the scaffold and genome alignment (+ means that the scaffold is spread among the genome)
    protected boolean missStart = false; // If this alignment shows a missassembly at the ContigStart attribute, it is set to true (Note: must NOT use reverse coord flag!)
    protected boolean missEnd = false; // If this alignment shows a missassembly at the ContigEnd attribute, it is set to true (Note: must NOT use reverse coord flag!)
    protected boolean genStart = false; // If this alignment is genuine at the ContigStart attribute, it is set to true (Note: must NOT use reverse coord flag!)
    protected boolean genEnd = false; // If this alignment is genuine at the ContigEnd attribute, it is set to true
    public boolean wasUsed = false; // if set to true, this scaffold was used to construct a super scaffold
    
    public AlignmentLine(String chr, int contigStart, int contigEnd, int contigSize, int alignStart, int alignEnd, String contigName){
        // This constructor is used by subclasses to set the needed attributes of the class without losing type casting abilities
        this.chr = chr;
        this.contigStart = contigStart;
        this.contigEnd = contigEnd;
        this.start = alignStart;
        this.end = alignEnd;
        this.name = contigName;
        this.contigSize = contigSize;
        
        this.percGenomeAlign = -1.0f;
        this.percSegmentAlign = -1.0f;
        this.proportionSegAlign = -1;
    }
    
    public AlignmentLine(String line){
        line = line.trim();
        String[] segs = line.split("\\|");
        
        //Genome Chr field
        segs[9] = segs[9].trim();
        Pattern chrnum = Pattern.compile("GK0{4,5}([0-9X]{1,2})\\.2");
        Matcher match = chrnum.matcher(segs[9]);
        if(match.find()){
            this.chr = "chr" + match.group(1).replaceFirst("^0", "");
        }else{
            try {
                throw new Exception("Chr match failed! " + segs[9]);
            } catch (Exception ex) {
                Logger.getLogger(AlignmentLine.class.getName()).log(Level.SEVERE, null, ex);
            }
        }
        
        //Contig name field
        segs[10] = segs[10].trim();
        this.name = segs[10];
        
        //Genome chr start and end
        String[] subsegs = segs[0].trim().split("\\s+");
        this.start = Integer.valueOf(subsegs[0]);
        this.end = Integer.valueOf(subsegs[1]);
        
        //Contig start and end
        subsegs = segs[1].trim().split("\\s+");
        this.contigStart = Integer.valueOf(subsegs[0]);
        this.contigEnd = Integer.valueOf(subsegs[1]);
        
        //Proportion seq align
        subsegs = segs[2].trim().split("\\s+");
        this.proportionSegAlign = Integer.valueOf(subsegs[1]) - Integer.valueOf(subsegs[0]);
        
        //Percentage of segment alignment
        segs[3] = segs[3].trim();
        this.percSegmentAlign = Float.valueOf(segs[3]);
        
        //Contig size information
        subsegs = segs[4].trim().split("\\s+");
        this.contigSize = Integer.valueOf(subsegs[1]);
        
        //Percentage of contig alignment to the genome
        subsegs = segs[5].trim().split("\\s+");
        this.percGenomeAlign = Float.valueOf(subsegs[1]);
        
               
    }
    
    public void setMissassembly(boolean start){
        if(start){
            this.missStart = true;
        }else{
            this.missEnd = true;
        }
    }
    
    public void setGenuine(boolean start){
        if(start)
            this.genStart = true;
        else
            this.genEnd = true;
    }
    
    /*
     * Getters
     */
    public int ContigStart(){
        return this.contigStart;
    }
    public int ContigEnd(){
        return this.contigEnd;
    }
    public float PercGenomeAlign(){
        return this.percGenomeAlign;
    }
    public float PercSegmentAlign(){
        return this.percSegmentAlign;
    }
    public int ProportionSeqAlign(){
        return this.proportionSegAlign;
    }
    public int ContigSize(){
        return this.contigSize;
    }
    
    public boolean equals(AlignmentLine t){
        return this.start == t.Start() && this.end == t.End() && this.Chr().equals(t.Chr()) && this.contigStart == t.contigStart && this.contigEnd == t.contigEnd;
    }
    
    @Override
    public int compareTo(BedAbstract o) {
        
            return this.Start() - o.Start();
        
    }
    
    
    /*
     * Testers
     */
    public boolean reversedCoords(){
        if(this.contigStart - this.contigEnd > 0){
            return true;
        }
        return false;
    }
    public boolean isMissassembly(){
        return (this.missEnd || this.missStart)? true : false;
    }
    public boolean MissStart(){
        return this.missStart;
    }
    public boolean MissEnd(){
        return this.missEnd;
    }
    
    /*
     * Printers
     */
    
    @Override
    public String toString(){
        StringBuilder str = new StringBuilder();
        str.append(this.chr).append("\t").append(this.start).append("\t");
        str.append(this.end).append("\t").append(this.name).append("\t");
        str.append(this.contigStart).append("\t").append(this.contigEnd).append("\t");
        str.append(this.proportionSegAlign).append("\t").append(this.percGenomeAlign);
        str.append("\t").append(this.percSegmentAlign);
        
        return str.toString();
    }
}

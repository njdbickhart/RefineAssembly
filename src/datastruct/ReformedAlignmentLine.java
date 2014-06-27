/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package datastruct;

/**
 *
 * @author bickhart
 */
public class ReformedAlignmentLine extends AlignmentLine{

    
    public ReformedAlignmentLine(String chr, int contigStart, int contigEnd, int contigSize, int alignStart, int alignEnd, String contigName){
        super(chr, contigStart, contigEnd, contigSize, alignStart, alignEnd, contigName);
    }
    
    public ReformedAlignmentLine(AlignmentLine a){
        super(a.Chr(), a.ContigStart(), a.ContigEnd(), a.ContigSize(), a.Start(), a.End(), a.Name());
    }
    
    @Override
    public String toString(){
        StringBuilder str = new StringBuilder();
        str.append(this.chr).append("\t").append(this.start).append("\t");
        str.append(this.end).append("\t").append(this.name).append("\t");
        str.append(this.contigStart).append("\t").append(this.contigEnd).append("\t");
        str.append(this.proportionSegAlign).append("\t").append(this.percGenomeAlign);
        str.append("\t").append(this.percSegmentAlign).append("\t").append("Reformed");
        
        return str.toString();
    }
}

/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package outputformats;

import datastruct.AGPChrMap;
import datastruct.AlignmentLine;
import datastruct.ContigFasta;
import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.charset.Charset;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.HashMap;

/**
 *
 * @author Dbick
 */
public class AGPOutput {
    private AGPChrMap data;
    
    public AGPOutput(AGPChrMap data){
        this.data = data;
    }
    
    public void generateAGP(String output, AGPFALogger log){
        try(BufferedWriter write = Files.newBufferedWriter(Paths.get(output), Charset.defaultCharset())){
            // Go through the chrs in karyotypic order
            for(String chr : utils.SortByChr.ascendingChr(data.getListChrs())){
                System.out.println("Generating AGP output for chr: " + chr);
                int chrPos = 1, ovCorr = 0;
                
                ArrayList<AlignmentLine> lines = this.data.getSortedAlignments(chr);
                
                for(int x = 0; x < lines.size(); x++){
                    AlignmentLine line = lines.get(x);
                    int astart = line.ContigStart(), aend = line.ContigEnd();

                    // Check if this alignment overlaps (or bookended) the last scaffold
                    // If so, modify the start coordinate to stitch the scaffolds together
                    // Otherwise, modify the start/end coordinates of the scaffold so that the full scaffold
                    // is incorporated into the AGP if there were no missassemblies on its ends.
                    if(ovCorr <= 0){
                        if(line.reversedCoords()){
                            astart = endCorrection(astart, line.MissStart(), log, line.ContigSize());
                            aend = startCorrection(aend, line.MissEnd(), log);
                        }else{
                            astart = startCorrection(astart, line.MissStart(), log);
                            aend = endCorrection(aend, line.MissEnd(), log, line.ContigSize());
                        }                        
                    }else{
                        if(line.reversedCoords()){
                            aend -= ovCorr;
                            astart = endCorrection(astart, line.MissStart(), log, line.ContigSize());
                        }else{
                            astart -= ovCorr;
                            aend = endCorrection(aend, line.MissEnd(), log, line.ContigSize());
                        }
                    }
                    
                    ovCorr = 0;
                    
                    int len = Math.abs(astart - aend);
                    String orient = (astart - aend > 0)? "-" : "+";

                    write.write(line.Chr() + "\t" + chrPos + "\t" + (chrPos + len) + "\t1\tP\t" +
                            line.Name() + "\t0\t" + astart + "\tscaffold\t" + aend + "\tno\t" +
                            orient + "\tna" + System.lineSeparator());
                    
                    chrPos += len;
                    
                    if(x < lines.size() -1 ){
                        // Check for a bookended or overlapping scaffold. If True, do not add an extra line for gaps
                        int ov = ovCount(line.Start(), line.End(), lines.get(x+1).Start(), lines.get(x+1).End());
                        if( ov > -1){
                            log.agpSummary("Scaffold overlapped next scaffold by " + ov + " bases : " + line.toString());
                            ovCorr = ov;
                        }else{
                            AlignmentLine next = lines.get(x+1);
                            int diff = next.Start() - chrPos;
                            
                            // Just to ensure that my alignments 
                            if(chrPos < next.Start()){
                                // I am only interested in big gaps of aligned sequence
                                if(diff > 20000){
                                    log.agpSummary("Chromosome alignment discontinuity. Previous alignment ended at " + chrPos + " and this alignment started at " + line.Start());
                                }
                                
                                // Since the next scaffold may extend beyond the alignment coordinates, check and reassign the diff value
                                int diffCorr = 0;
                                if(next.reversedCoords()){
                                    diffCorr = startCorrection(next.ContigEnd(), next.MissEnd(), log);
                                    diffCorr = next.ContigEnd() - diffCorr;
                                }else{
                                    diffCorr = startCorrection(next.ContigStart(), next.MissStart(), log);
                                    diffCorr = next.ContigStart() - diffCorr;
                                }
                                
                                // This should eliminate the chance that the difference ends up being negative for this region
                                // Since such regions were not represented in the alignment, I assume that they are actual differences
                                if(diff - diffCorr > 0){
                                    diff -= diffCorr;
                                }
                                
                                write.write(line.Chr() + "\t" + chrPos + "\t" + (chrPos + diff) + "\t1\tN\t" +
                                    "gap" + "\t" + diff + "\t" + chrPos + "\tscaffold\t" + (chrPos + diff) + "\tyes\t" +
                                    "+" + "\talign_genus" + System.lineSeparator());
                                chrPos = next.Start();
                            }else if(diff < 0){
                                write.write(line.Chr() + "\t" + chrPos + "\t" + (chrPos + 100) + "\t1\tU\t" +
                                    "gap" + "\t" + 100 + "\t" + chrPos + "\tscaffold\t" + (chrPos + 100) + "\tyes\t" +
                                    "+" + "\talign_genus" + System.lineSeparator());
                                chrPos += 100;
                                log.errorMessage("Next scaffold aligment start less than chrPos! " + chrPos + " " + next.toString());
                            }
                        }
                    }
                }
            }
        }catch(IOException ex){
            log.errorMessage("Problem writing to output: " + output);
            log.errorMessage(ex.getMessage());
        }
    }
    
    private int startCorrection(int astart, boolean startMiss, AGPFALogger log){
        int value;
        
        // If the start end of the contig does not align to the genome, and there are no misassemblies, then might as well include it!
        value = astart;
        if(astart < 10000 && !startMiss){
            log.agpSummary("Extending scaffold edge to include the beginning: " + astart );
            value = 1;
        }
       
        return value;
    }
    
    private int endCorrection(int aend, boolean endMiss, AGPFALogger log, int scafSize){
        int value;
        
        // If the end of the contig does not align to the genome, and there are no misassemblies, then might as well include it!
        value = aend;
        if(aend > scafSize - 10000 && !endMiss){
            log.agpSummary("Extending scaffold edge to include the end: " + aend + " " + scafSize);
            value = scafSize;
        }
     
        return value;
    }
    
    protected int ovCount (int start1, int end1, int start2, int end2){
        return soonest(end1, end2) - latest(start1, start2);
    }
    protected int soonest (int a, int b){
        return (a >= b)? b : a;
    }
    protected int latest (int a, int b){
        return (a >= b)? a : b;
    }
}

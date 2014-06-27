/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package datastruct;

import java.util.ArrayList;

/**
 *
 * @author Dbick
 */
public class ContigFasta {
    private String scaffold;
    private int linecount;
    private String sequence;
    
    public ContigFasta(String sequence, ArrayList<AlignmentLine> lines){
        this.sequence = sequence;
        this.scaffold = lines.get(0).Name();
        this.linecount = lines.size();
    }
    
    public String getSequenceFromScaffold(AlignmentLine line){
        if(!line.Name().equals(this.scaffold)){
            System.out.println("ERROR! AlignmentLine input does not match scaffold!");
            System.exit(-1);
        }
        
        int astart = line.contigStart, aend = line.contigEnd;
        
        // If the tail end of the contig does not align to the genome, and there are no misassemblies, then might as well include it!
        if(astart > line.contigSize - 20000 && !line.MissStart()){
            astart = line.contigSize;
        }else if(aend > line.contigSize - 20000 && !line.MissEnd()){
            aend = line.contigSize;
        }
        
        linecount--;
        
        if(astart > aend){
            return revComp(this.sequence.substring(aend, astart));
        }else{
            return this.sequence.substring(astart, aend);
        }
        
        
    }
    
    private String revComp(String subseq){
        char[] letters = subseq.toCharArray();
        char[] comp = new char[letters.length];
        
        for(int x = 0; x < letters.length; x++){
            int rev = letters.length - 1 - x;
            switch(letters[x]){
                case 'A':
                case 'a':
                    comp[rev] = 'T';
                    break;
                case 'G':
                case 'g':
                    comp[rev] = 'C';
                    break;
                case 'C':
                case 'c':
                    comp[rev] = 'G';
                    break;
                case 'T':
                case 't':
                    comp[rev] = 'A';
                    break;
                default:
                    comp[rev] = 'N';
                    break;
            }
        }
        return new String(comp);
    }
}

/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package factories;

import datastruct.AlignmentLine;
import java.io.BufferedReader;
import java.io.IOException;
import java.nio.charset.Charset;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Set;
import outputformats.alignmentLogger;

/**
 *
 * @author Dbick
 */
public class AlignmentLineFactory {
    private HashMap<String, Integer> contigList = new HashMap<>();
    
    //HashMap<Chromosome, Array of alignment lines>
    private HashMap<String, ArrayList<AlignmentLine>> contigGrouping = new HashMap<>();
    
    public void generateAlignmentData(String file, alignmentLogger log){
        try(BufferedReader input = Files.newBufferedReader(Paths.get(file), Charset.defaultCharset())){
            String line;
            while((line = input.readLine())!= null){
                AlignmentLine align = new AlignmentLine(line);
                updateContigList(line);
                if(!this.contigGrouping.containsKey(align.Chr())){
                    this.contigGrouping.put(align.Chr(), new ArrayList());
                }
                this.contigGrouping.get(align.Chr()).add(align);
            }
        }catch(IOException ex){
            log.errorMessage("Input error in Alignment Line Factory!");
        }
    }
    
    private void updateContigList(String line){
        line = line.trim();
        String[] segs = line.split("\\|");
        String key = segs[segs.length - 1].trim();
        int i = 1;
        if(contigList.containsKey(key)){
            i += contigList.get(key);
        }
        contigList.put(key, i);
    }
    
    /*
     * Getters
     */
    public Set<String> getContigList(){
        return this.contigList.keySet();
    }
    public int getContigCount(String contigname){
        return this.contigList.get(contigname);
    }
    public Set<String> getChrList(){
        return this.contigGrouping.keySet();
    }
    public ArrayList<AlignmentLine> getAlignFromChr(String chr){
        return this.contigGrouping.get(chr);
    }
    public ArrayList<AlignmentLine> collectAlignFromContig(String contigname){
        ArrayList<AlignmentLine> collection = new ArrayList<>();
        for(String chr : this.contigGrouping.keySet()){
            for(AlignmentLine a : this.contigGrouping.get(chr)){
                if(a.Name().equals(contigname)){
                    collection.add(a);
                }
            }
        }
        return collection;
    }
    public HashMap<String, ArrayList<AlignmentLine>> getAllAligns(){
        return this.contigGrouping;
    }
}

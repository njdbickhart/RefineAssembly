/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package workers;

import datastruct.AlignmentLine;
import datastruct.ReformedAlignmentLine;
import factories.AlignmentLineFactory;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.charset.Charset;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Set;
import java.util.TreeSet;
import outputformats.alignmentLogger;
import utils.LineIntersect;

/**
 *
 * @author Dbick
 */
public class mateWindowConsistency {
    /*****************************/
    // In order to refine the distance from which a zero-mate window is detected, modify these two variables:
    /*****************************/
    private final int zeroWinIncrement = 500; // The size increment for the zero window detection algorithm
    private final int maxSearchCount = 1; // The counter of times to increment the search distance in the zero window detection algorithm
    protected HashMap<String, ArrayList<AlignmentLine>> contigGrouping = new HashMap<>(); // From AlignmentLineFactory; {scaffold}->array
    protected ArrayList<ReformedAlignmentLine> finished = new ArrayList<>(); // Container for fixed alignments
    protected String mateFile;
    protected ArrayList<String> genuine = new ArrayList<>();
    protected HashSet<String> missassemblies = new HashSet<>();
    protected ArrayList<Integer> distances = new ArrayList<>(); // Container for distances of zero-mate alignments from alignment breakpoints
    
    public mateWindowConsistency(AlignmentLineFactory contigGrouping, alignmentComparison comparison, String mateFile){
        Set<String> problems = comparison.getAllProblemScaffolds();
        
        for(String scaffold : problems){
            ArrayList<AlignmentLine> temp = contigGrouping.collectAlignFromContig(scaffold);
            this.contigGrouping.put(scaffold, temp);
        }
        
        this.mateFile = mateFile;
    }
    
    public void searchMatesResolveErrors(alignmentLogger log){
        String line = null;
        try(BufferedReader reader = Files.newBufferedReader(Paths.get(this.mateFile), Charset.defaultCharset())){
            // Just to avoid the penalties of creating new classes every time, I will use an old C trick here with array indicies
            // Here are the arrays for the mate pair mappings
            //ArrayList<Integer> positions = new ArrayList<>();
            //ArrayList<Integer> matecounts = new ArrayList<>();
            
            //findAllZeroMateLocs(positions, matecounts, progress, total, log, corrected, reader);
            findEdgeZeroMateLocs(reader, log);
            
        }catch(IOException ex){
            log.errorMessage("Error Reading Mate pair file: " + this.mateFile + "! Terminating program prematurely...");
        }catch(java.lang.ArrayIndexOutOfBoundsException ex){
            System.err.println("Error at line: " + line);
            ex.printStackTrace();
        }
        
    }
    
    public ArrayList<String> getGenuineScaffolds(){
        return this.genuine;
    }
    
    protected synchronized void addToFinal(ReformedAlignmentLine a){
        this.finished.add(a);
    }
    
    private ArrayList<AlignmentLine> getOverlappingEntries(ArrayList<AlignmentLine> lines, int start, int end){
        ArrayList<AlignmentLine> holder = new ArrayList<>();
        for(AlignmentLine a : lines){
            if(LineIntersect.ovCount(a.ContigStart(), a.ContigEnd(), start, end) > 1){
                holder.add(a);
            }
        }
        return holder;
    }
    
    private float Average(ArrayList<Integer> list){
        float sum = 0;
        if(list.isEmpty()){return 0.0f;}
        for(int a : list){
            sum += a;
        }
        return sum / list.size();
    }
    
    private double StDev(ArrayList<Integer> list, float avg){
        float stdev = 0.0f;
        if(list.isEmpty()){return 0.0f;}
        
        for(int x = 0; x < list.size(); x++){
            stdev += (float) Math.pow(list.get(x) - avg, 2.0f);
        }
        double variance = stdev / (float) (list.size() - 1);
        return Math.sqrt(variance);
    }
    
    private float avgStep(ArrayList<Integer> list){
        ArrayList<Integer> values = new ArrayList<>(list.size());
        for(int x = 0; x < list.size() - 1; x++){
            values.add(Math.abs(list.get(x+1) - list.get(x)));
        }
        return Average(values);
    }

    private void findEdgeZeroMateLocs(BufferedReader reader, alignmentLogger log) throws IOException{
        // This method identifies the edges of problematic alignment lines and searches for mate pair concordancy
        // If there are zero mate pair alignments nearby, then the alignment is from a misassembly
        String line, scaffold = null;
        int progress = 0, total = this.contigGrouping.keySet().size(), uncorrected = 0, corrected = 0;
        ArrayList<Integer> positions = new ArrayList<>();
        ArrayList<Integer> matecounts = new ArrayList<>();
        
        while((line = reader.readLine()) != null){
            String[] segs = line.split("\\s+");
            if(scaffold == null){
                scaffold = segs[0];
            }
            if(segs.length < 3){continue;}
            
            // This is my really simple way of getting all of the data and working on the results
            // It assumes that the mate pair alignment file is sorted with respect to the contig names
            if(!scaffold.equals(segs[0]) && this.contigGrouping.containsKey(scaffold)){
                progress++;
                System.out.println("Working on contig: " + progress + " of total: " + total);
                ArrayList<AlignmentLine> working = this.contigGrouping.get(scaffold);
                
                //log.alignmentSummary("mateWindowConsistency", "Read Positions and mate counts for scaffold: " + scaffold + ".");
                
                // Calculate the average window size for the matepair position file
                float avgWinSize = 700;
                
                // Only work on the scaffold if it has multiple alignments
                if(working.size() > 1){
                    Collections.sort(working);                    
                    for(AlignmentLine a : working){
                        // I need to check both sides of the alignment line if the alignment is not near the scaffold start or end
                        // if the alignment start is near the scaffold start then skip that side
                        if(!(a.reversedCoords() && a.ContigEnd() < avgWinSize) || !(a.ContigStart() < avgWinSize)){
                            if(a.reversedCoords()){
                                boolean fix = checkForZeroMates(positions, matecounts, a.ContigEnd(), avgWinSize);
                                if(fix){
                                    a.setMissassembly(false);
                                    this.missassemblies.add(a.Name());
                                    log.alignmentSummary("mateWindowConsistency1", "Missassembly! Scaffold: " + scaffold + " end: " + a.ContigEnd());
                                }else{
                                    a.setGenuine(false);
                                    log.alignmentSummary("mateWindowConsistency1", "Genuine difference link found: Scaffold: " + scaffold + " end: " + a.ContigEnd());
                                }
                            }else{
                                boolean fix = checkForZeroMates(positions, matecounts, a.ContigStart(), avgWinSize);
                                if(fix){
                                    a.setMissassembly(true);
                                    this.missassemblies.add(a.Name());
                                    log.alignmentSummary("mateWindowConsistency2", "Missassembly! Scaffold: " + scaffold + " start: " + a.ContigStart());
                                }else{
                                    a.setGenuine(true);
                                    log.alignmentSummary("mateWindowConsistency2", "Genuine difference link found: Scaffold: " + scaffold + " start: " + a.ContigStart());
                                }
                            }
                        }
                        
                        // if the alignment end is near the scaffold end, then skip that side
                        if((a.reversedCoords() && a.ContigStart() < (a.ContigSize() - avgWinSize)) || (!a.reversedCoords() && a.ContigEnd() < (a.ContigSize() - avgWinSize))){
                            if(a.reversedCoords()){
                                boolean fix = checkForZeroMates(positions, matecounts, a.ContigStart(), avgWinSize);
                                if(fix){
                                    a.setMissassembly(true);
                                    this.missassemblies.add(a.Name());
                                    log.alignmentSummary("mateWindowConsistency3", "Missassembly! Scaffold: " + scaffold + " start: " + a.ContigStart());
                                }else{
                                    a.setGenuine(true);
                                    log.alignmentSummary("mateWindowConsistency3", "Genuine difference link found: Scaffold: " + scaffold + " start: " + a.ContigStart());
                                }
                            }else{
                                boolean fix = checkForZeroMates(positions, matecounts, a.ContigEnd(), avgWinSize);
                                if(fix){
                                    a.setMissassembly(false);
                                    this.missassemblies.add(a.Name());
                                    log.alignmentSummary("mateWindowConsistency4", "Missassembly! Scaffold: " + scaffold + " end: " + a.ContigEnd());
                                }else{
                                    a.setGenuine(true);
                                    log.alignmentSummary("mateWindowConsistency4", "Genuine difference link found: Scaffold: " + scaffold + " end: " + a.ContigEnd());
                                }
                            }
                        }
                        
                        if(!a.isMissassembly()){
                            // If we did not detect a missassembly, then keep track of the number of entries that were actual differences
                            log.alignmentSummary("mateWindowConsistency5", "Full length, genuine scaffold alignment candidate: " + scaffold + " " + a.toString());
                        }
                    }
                    boolean isMiss = false;
                    for(AlignmentLine a : working){
                        if( a.isMissassembly()){
                            isMiss = true;
                        }
                    }
                    if(isMiss){
                        corrected++;
                    }else{
                        this.genuine.add(scaffold);
                        uncorrected++;                        
                    }
                }
                
                positions.clear();
                matecounts.clear();
                scaffold = segs[0];
            }else if(!scaffold.equals(segs[0]) && !this.contigGrouping.containsKey(scaffold)){
                positions.clear();
                matecounts.clear();
                scaffold = segs[0];
            } 
            positions.add(Integer.valueOf(segs[1]));
            matecounts.add(Integer.valueOf(segs[2]));
        }
        log.alignmentSummary("mateWindowConsistency", "Identified: " + uncorrected + " genuine differences and " + corrected + " missassemblies.");
    }
    
    public void printOutGenuineChromosomeDiffs(String outfile){
        try(BufferedWriter output = Files.newBufferedWriter(Paths.get(outfile), Charset.defaultCharset())){
            for(String scaffold : this.genuine){
                for(AlignmentLine l : this.contigGrouping.get(scaffold)){
                    output.write(scaffold + "\t" + l.toString() + System.lineSeparator());
                }
            }
        }catch(IOException ex){
            ex.printStackTrace();
        }
    }
    
    public void printOutMissassemblies(String outfile){
        try(BufferedWriter output = Files.newBufferedWriter(Paths.get(outfile), Charset.defaultCharset())){
            for(String scaffold : this.missassemblies){
                for(AlignmentLine l : this.contigGrouping.get(scaffold)){
                    output.write(scaffold + "\t" + l.toString() + System.lineSeparator());
                }
            }
        }catch(IOException ex){
            ex.printStackTrace();
        }
    }
    
    private boolean checkForZeroMates(ArrayList<Integer> positions, ArrayList<Integer> matecounts, int position, float avgWinSize){
        // Start at the second alignment and stop before the end; the start and end tend to have fewer mate pair mappings by default
        int beginning = position - zeroWinIncrement;
        int ending = position + zeroWinIncrement;
        int StopCounter = 0;
        // The following array is a list of array indicies that have low mate coverage
        ArrayList<Integer> intersecting = new ArrayList<>();
        
        while(beginning > 0 && ending < positions.get(positions.size() -1)){
            if(StopCounter >= maxSearchCount) break;
            
            // internal loop to search through the beginning/end coordinates
            for(int x = 1; x < positions.size() - 1; x++){
                // Find the window indicies that contain the event
                if(positions.get(x) >= beginning && positions.get(x-1) < ending){
                    // Now that we've found the position within the windows, check to see if we have a breakpoint
                    // Criteria are to find a "slope" of low coverage regions
                    if(matecounts.get(x) <= 1 && matecounts.get(x -1) <= 3 && matecounts.get(x+1) <= 3){
                        intersecting.add(x);
                    }else if(matecounts.get(x) <= 4 && matecounts.get(x - 1) <= 4 && positions.get(x) - positions.get(x - 1) > avgWinSize){
                        intersecting.add(x);
                    }
                }
            }
            
            // If we have intersections, loop through them to find the closest one
            if(!intersecting.isEmpty()){
                int size = 1000000;
                int index = 0;
                for(int i : intersecting){
                    if(size > Math.abs(position - positions.get(i))){
                        index = i;
                        int current = positions.get(i);
                        int value = matecounts.get(i);
                        size = Math.abs(position - current);
                    }                    
                }
                this.distances.add(positions.get(index) - position);
                return true;
            }                
            
            // Increment search coordinates
            beginning -= zeroWinIncrement;
            ending += zeroWinIncrement;
            StopCounter++;
        }
        return false;
    }

    @Deprecated
    private void findAllZeroMateLocs(ArrayList<Integer> positions, ArrayList<Integer> matecounts, int progress, int total, alignmentLogger log, int corrected, BufferedReader reader) throws NumberFormatException, IOException {
        String line;
        // Populate the mate pair arrays and calculate statistics
        String scaffold = null;
        while((line = reader.readLine()) != null){
            String[] segs = line.split("\\s+");
            if(scaffold == null){
                scaffold = segs[0];
            }
            if(segs.length < 3){continue;}
            
            positions.add(Integer.valueOf(segs[1]));
            matecounts.add(Integer.valueOf(segs[2]));
            
            if(!scaffold.equals(segs[0]) && this.contigGrouping.containsKey(scaffold)){
                
                progress++;
                System.out.println("Working on contig: " + progress + " of total: " + total);
                ArrayList<AlignmentLine> working = this.contigGrouping.get(scaffold);
                
                log.alignmentSummary("mateWindowConsistency", "Read Positions and mate counts for scaffold: " + scaffold + ".");

                float avg = Average(matecounts);
                float stdev = (float)StDev(matecounts, avg);


                // Calculate the average window size for the matepair position file
                float avgWinSize = Average(positions);


                // Loop through the mate arrays looking for mate pair mapping inconsistencies
                int lowerthresh = 1;
                boolean appliedCorrection = false;
                int lastposition = 1;
                for(int x = 1; x < matecounts.size(); x++){

                    if(matecounts.get(x) <= lowerthresh){
                        // Look ahead to see if there is a lower count of mate pairs ahead                    
                        int z = x; // counter to hold last position of y
                        for(int y = x + 1; y < matecounts.size(); y++){
                            if(matecounts.get(y) > matecounts.get(x)){
                                break;
                            }else if(matecounts.get(y) < matecounts.get(x) && matecounts.get(y) < matecounts.get(z)){
                                z = y;
                            }
                        }

                        if(z != x){
                            // We found a region that had lower matecounts than x
                            x = z;
                        }

                        // Now find the alignment lines that this influences
                        ArrayList<AlignmentLine> influenced = getOverlappingEntries(working, lastposition, positions.get(x));
                        for(AlignmentLine a : influenced){
                            if((a.reversedCoords() && a.ContigStart() > positions.get(x))){
                                corrected++;
                                this.finished.add(new ReformedAlignmentLine(a.Chr(), positions.get(x), a.ContigEnd(), a.ContigSize(), a.Start(), a.End(), a.Name()));
                                log.alignmentSummary("mateWindowConsistency", "Scaffold: " + a.Name() + " missassembly. Setting start to: " + positions.get(x) + " from: " + a.ContigStart() + ".");
                            }else if(!a.reversedCoords() && a.ContigEnd() > positions.get(x)){
                                this.finished.add(new ReformedAlignmentLine(a.Chr(), a.ContigStart(), positions.get(x), a.ContigSize(), a.Start(), a.End(), a.Name()));
                                log.alignmentSummary("mateWindowConsistency", "Scaffold: " + a.Name() + " missassembly. Setting end to: " + positions.get(x) + " from: " + a.ContigEnd() + ".");
                            }else{
                                this.finished.add(new ReformedAlignmentLine(a));
                            }
                        }
                        
                        // Look ahead again to make sure that we aren't doing a separate analysis for contiguous low mate regions
                        z = x;
                        for(int y = x + 1; y < matecounts.size(); y++){
                            if(matecounts.get(y) > lowerthresh){
                                break;
                            }else if(matecounts.get(y) <= matecounts.get(x) && matecounts.get(y) < matecounts.get(z)){
                                z = y;
                            }
                        }
                        
                        if(z != x){
                            // We found a region that had lower matecounts than x
                            x = z;
                        }
                        
                        lastposition = positions.get(x);
                    }               

                }

                if(!appliedCorrection){
                    // The mate pair counts were fine for this scaffold, so we can put its alignment lines into the good bin
                    for(AlignmentLine a : working){
                        this.finished.add(new ReformedAlignmentLine(a));
                    }
                }
                positions.clear();
                matecounts.clear();
                scaffold = segs[0];
            }else if(!scaffold.equals(segs[0]) && !this.contigGrouping.containsKey(scaffold)){
                scaffold = segs[0];
            }
        }
    }
    
    /*
     * Getters
     */
    public ArrayList<String> getGenuine(){
        return this.genuine;
    }
    public HashSet<String> getMissassemblies(){
        return this.missassemblies;
    }
    public ArrayList<Integer> getBreakPoints(){
        return this.distances;
    }
}

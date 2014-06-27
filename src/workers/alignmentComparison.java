/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package workers;

import datastruct.AGPChrMap;
import datastruct.AlignmentLine;
import factories.AlignmentLineFactory;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.Set;
import java.util.TreeSet;
import outputformats.alignmentLogger;

/**
 *
 * @author Dbick
 */
public class alignmentComparison {
    private HashSet<String> crossChrProblemContigs = new HashSet<>();
    private HashSet<String> intraScaffoldProblemContigs = new HashSet<>();
    private AGPChrMap normalScaffolds = new AGPChrMap();
    private alignmentLogger log;
    
    public alignmentComparison(alignmentLogger log){
        this.log = log;
    }
    
    
    public void crossChromosomeComparison(AlignmentLineFactory factory){
        int problems = 0;
        for(String contig : factory.getContigList()){
            ArrayList<AlignmentLine> alignments = factory.collectAlignFromContig(contig);
            Set<String> chrSet = new TreeSet<>();
            for(AlignmentLine a : alignments){
                chrSet.add(a.Chr());
            }
            if(chrSet.size() > 1){
                crossChrProblemContigs.add(contig);
                problems++;
            }
        }
        
        log.alignmentSummary("alignmentComparison", "Identified " + problems + " cross chromosome problems out of " + factory.getContigList().size() + " total contigs");
    }
    
    public void intraScaffoldComparison(AlignmentLineFactory factory){
        int problems = 0;
        for(String contig : factory.getContigList()){
            if(this.crossChrProblemContigs.contains(contig)){
                // I need to deal with this contig a different way
                continue;
            }
            ArrayList<AlignmentLine> alignments = factory.collectAlignFromContig(contig);
            int distance = 0;
            if(alignments.size() > 1){
                int lastend = 0;
                for(AlignmentLine a : alignments){
                    if(lastend == 0){
                        if(a.reversedCoords()){
                            lastend = a.Start();
                        }else{
                            lastend = a.End();
                        }
                        continue;
                    }
                    if(a.reversedCoords()){
                        distance += (lastend > a.End())? lastend - a.End() : a.End() - lastend;
                    }else{
                        distance += (lastend > a.Start())? lastend - a.Start() : a.Start() - lastend;
                    }
                }
            }
            if(distance > 1000000){
                this.intraScaffoldProblemContigs.add(contig);
                problems++;
            }
            
        }
        log.alignmentSummary("alignmentComparison", "Identified " + problems + " intrascaffold problems out of " + factory.getContigList().size() + " total contigs");
    }
    
    public void arrangeNormalScaffolds(AlignmentLineFactory factory){
        int count = 0;
        for(String contig : factory.getContigList()){
            if(this.crossChrProblemContigs.contains(contig) || this.intraScaffoldProblemContigs.contains(contig)){
                continue;
            }else{
                ArrayList<AlignmentLine> alignments = factory.collectAlignFromContig(contig);
                for(AlignmentLine a : alignments){
                    this.normalScaffolds.addBedData(a);
                    count++;
                }
            }
        }
        log.alignmentSummary("alignmentComparison", "Placed " + count + " scaffolds into a normal holding container.");
    }
    
    public void printChrConflictList(AlignmentLineFactory factory){
        for(String contig : crossChrProblemContigs){
            System.out.println("Chr Conflict: " + contig);
            ArrayList<AlignmentLine> problems = factory.collectAlignFromContig(contig);
            for(AlignmentLine a : problems){
                System.out.println(a);
            }
        }
    }
    
    public void printIntraConflictList(AlignmentLineFactory factory){
        for(String contig: intraScaffoldProblemContigs){
            System.out.println("Scaff Conflict: " + contig);
            ArrayList<AlignmentLine> problems = factory.collectAlignFromContig(contig);
            for(AlignmentLine a : problems){
                System.out.println(a);
            }
        }
    }
    
    
    /*
     * Getters
     */
    
    public Set<String> getAllProblemScaffolds(){
        Set<String> problems = new HashSet<>();
        problems.addAll(this.crossChrProblemContigs);
        problems.addAll(this.intraScaffoldProblemContigs);
        return problems;
    }
    
    public Set<String> getCrossChrProblemScaffolds(){
        return this.crossChrProblemContigs;
    }
    
    public Set<String> getIntraScaffoldProblems(){
        return this.intraScaffoldProblemContigs;
    }
}

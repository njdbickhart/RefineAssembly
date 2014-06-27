/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package refineassemblydraft;

import datastruct.AGPChrMap;
import factories.AlignmentLineFactory;
import factories.SuperScaffoldBuilder;
import outputformats.AGPFALogger;
import outputformats.AGPOutput;
import outputformats.DistanceZeroMateOut;
import outputformats.alignmentLogger;
import workers.alignmentComparison;
import workers.mateWindowConsistency;

/**
 *
 * @author Dbick
 */
public class RefineAssemblyDraft {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // Process command line arguments and create log file
        ParseCommandLine cmd = new ParseCommandLine(args);
        alignmentLogger log = new alignmentLogger(cmd.alignmentLog);
        
        // Read input nucmer alignment lines
        AlignmentLineFactory alignFac = new AlignmentLineFactory();
        alignFac.generateAlignmentData(cmd.alignmentFile, log);
        
        // Identify putative misassemblies from alignments
        alignmentComparison comparisonTool = new alignmentComparison(log);
        comparisonTool.crossChromosomeComparison(alignFac);
        comparisonTool.intraScaffoldComparison(alignFac);
        
        // Debug: alignment test
        //comparisonTool.printChrConflictList(alignFac);
        //comparisonTool.printIntraConflictList(alignFac);
        
        // Read information from the mate pair file to confirm misassemblies
        mateWindowConsistency mateCheck = new mateWindowConsistency(alignFac, comparisonTool, cmd.matepairFile);
        mateCheck.searchMatesResolveErrors(log);
        mateCheck.printOutGenuineChromosomeDiffs(cmd.output + ".diffs");
        mateCheck.printOutMissassemblies(cmd.output + ".miss");
        
        // Diagnostics: print out the distances table
        DistanceZeroMateOut mateOut = new DistanceZeroMateOut();
        mateOut.PrintOutData(mateCheck.getBreakPoints(), cmd.output + ".dist.tab", 100);
        
        // Generate contiguous chromosome fastas and AGP file containing scaffold
        //AGPFALogger agpLog = new AGPFALogger(cmd.agpFaLog);
        //AGPChrMap alignMap = new AGPChrMap();
        //alignMap.addAlignmentLines(alignFac);
        
        //AGPOutput agp = new AGPOutput(alignMap);
        //agp.generateAGP(cmd.output + ".agp", agpLog);
        
        // Create super scaffolds and print them out
        SuperScaffoldBuilder ssBuild = new SuperScaffoldBuilder(alignFac, mateCheck);
        ssBuild.mergeAlignments();
        ssBuild.seedExtendScaffolds();
        //ssBuild.chainScaffolds();
        //ssBuild.linkUpScaffolds();
        ssBuild.printOutScaffoldPlan(cmd.output + ".sscafs");
        
        // Close files and exit the program
        log.close();
        //agpLog.close();
    }
}

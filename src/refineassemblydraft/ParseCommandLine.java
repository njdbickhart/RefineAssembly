/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package refineassemblydraft;

/**
 *
 * @author bickhart
 */
public class ParseCommandLine {
    public String nl = System.lineSeparator();
    public String usage = "USAGE: java -jar RefineAssemblyDraft.jar" + nl
            + "\t-a\tInput Nucmer alignment file" + nl
            + "\t-m\tMate pair mapping count file" + nl
            + "\t-o\tbase output name" + nl
            + nl
            + "[Optional arguments]" + nl
            + "\t-al\tAlignment log file name [alignment.log]"
            + "\t-fl\tAGP/FA log file name [agpFa.log]";
    public String alignmentFile;
    public String alignmentLog = "alignment.log";
    public String agpFaLog = "agpFa.log";
    public String matepairFile;
    public String output;
    
    
    public ParseCommandLine(String[] args){
        for(int i = 0; i < args.length; i++){
            switch(args[i]){
                case "-a":
                    this.alignmentFile = args[i+1];
                    break;
                case "-al":
                    this.alignmentLog = args[i+1];
                    break;
                case "-fl":
                    this.agpFaLog = args[i+1];
                    break;
                case "-m":
                    this.matepairFile = args[i+1];
                    break;
                case "-o":
                    this.output = args[i+1];
                    break;
            }
        }
    }
}

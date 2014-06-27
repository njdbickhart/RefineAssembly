/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package outputformats;

import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.charset.Charset;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.sql.Timestamp;
import java.util.Date;

/**
 *
 * @author bickhart
 */
public class AGPFALogger {
    private BufferedWriter logHandle;
    
    public AGPFALogger(String logName){
        try{
            logHandle = Files.newBufferedWriter(Paths.get(logName), Charset.defaultCharset());
        }catch(IOException ex){
            System.err.println("Could not create log file: " + logName + "! Exiting...");
        }
    }
    
    /*
     * output methods
     */
    public synchronized void errorMessage(String msg){
        try{
            Date date = new Date();
            this.logHandle.write("ERROR -" + new Timestamp(date.getTime()) + " - " + msg);
            this.logHandle.write(System.lineSeparator());
        }catch(IOException ex){
            System.err.println("Error writing to log file: " + this.logHandle.toString() + " for error message: " + msg);
        }
    }
    
    public synchronized void agpSummary(String msg){
        try{
            Date date = new Date();
            this.logHandle.write("AGP - " + new Timestamp(date.getTime()) + " - " + msg);
            this.logHandle.write(System.lineSeparator());
        }catch(IOException ex){
            System.err.println("Error writing");
        }
    }
    
    public synchronized void faSummary(String msg){
        try{
            Date date = new Date();
            this.logHandle.write("FA *** " + new Timestamp(date.getTime()) + " *** " + msg);
            this.logHandle.write(System.lineSeparator());
        }catch(IOException ex){
            System.err.println("Error writing");
        }
    }
    
    /*
     * File control methods
     */
    public void close(){
        try{
            this.logHandle.close();
        }catch(IOException ex){
            System.err.println("Could not close alignment log file: " + this.logHandle.toString() + "!");
        }
    }
}

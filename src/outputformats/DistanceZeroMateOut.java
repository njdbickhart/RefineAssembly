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
import java.util.ArrayList;
import java.util.Collections;

/**
 *
 * @author bickhart
 */
public class DistanceZeroMateOut {
    
    public void PrintOutData(ArrayList<Integer> data, String output, int binValue){
        Collections.sort(data);
        int min = data.get(0);
        int max = data.get(data.size() - 1);
        
        // rounding bin values to the outer ends of the min and max coordinates
        int remainder = min % binValue;
        min -= binValue - Math.abs(remainder);
        
        remainder = max %binValue;
        max += binValue - remainder;
        
        int divisor = (max - min) / binValue;
        int[] indicies = new int[divisor + 1];
        int[] values = new int[divisor + 1];
        
        int t = 0;
        for(int x = min; x <= max; x += binValue){
            indicies[t] = x;
            t++;
        }
        for(int i : data){
            for(int x = 0; x < indicies.length; x++){
                int start = indicies[x];
                int end = indicies[x] + binValue;
                
                if(i >= start && i < end){
                    values[x] += 1;
                }
            }
        }
        
        // Now to print it all out
        try(BufferedWriter out = Files.newBufferedWriter(Paths.get(output), Charset.defaultCharset())){
            for(int x = 0; x < indicies.length; x++){
                out.write(indicies[x] + "\t" + values[x] + "\n");
            }
        }catch(IOException ex){
            ex.printStackTrace();
        }
    }
}

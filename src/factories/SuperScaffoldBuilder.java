/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package factories;

import datastruct.AlignmentLine;
import datastruct.SuperScaffold;
import file.BedAbstract;
import file.BedMap;
import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.charset.Charset;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.TreeMap;
import workers.mateWindowConsistency;

/**
 *
 * @author desktop
 */
public class SuperScaffoldBuilder {
    private HashMap<String, ArrayList<AlignmentLine>> bArray; // hash that lists alignment lines by Chr
    private HashMap<String, ArrayList<AlignmentLine>> sArray; // hash that lists alignment lines by scaffold name
    private ArrayList<AlignBed> mergeStore;
    private ArrayList<SuperScaffold> scaffolds; 
    private HashSet<String> genuines;
    private HashSet<String> misassemblies;
    private ArrayList<Set<String>> alignStore = new ArrayList<>();
    private Set<String> mergedSample = new HashSet<String>();
    private HashMap<String, ArrayList<AlignBed>> alignTree = new HashMap<>();
    private int ovlpBaseNum = 50000;
    
    public SuperScaffoldBuilder(AlignmentLineFactory fact, mateWindowConsistency mates){
        this.sArray = new HashMap<>();
        this.bArray = fact.getAllAligns();
        this.alignLineHasher(fact);
        
        this.genuines = new HashSet<>();
        for(String scaf : mates.getGenuine()){
            this.genuines.add(scaf);
        }
        this.misassemblies = mates.getMissassemblies();
        this.scaffolds = new ArrayList<>();
    }
    
    // First step! Take all scaffolds that align to the genome within 50bp of each other (or overlapping) and group them
    public void mergeAlignments() {
        int merger = 0;
        this.mergeStore = new ArrayList<AlignBed>();
        for(String chr : this.bArray.keySet()){
            // Sorting
            ArrayList<AlignmentLine> working = this.bArray.get(chr);
            Collections.sort(working, new Comparator(){

                @Override
                public int compare(Object t, Object t1) {
                    AlignmentLine mine = (AlignmentLine) t;
                    AlignmentLine yours = (AlignmentLine) t1;
                    return mine.Start() - yours.Start();
                }
                
            });

            //Now the merger
            
            ArrayList<AlignmentLine> aStorage = new ArrayList<>();
            //aStorage.add(working.get(0));
            AlignmentLine prev;
            boolean overlap = false; // boolean counter to signify that overlapping is taking place
            int minStart = -1;
            int maxEnd = -1;
            
            for(int i = 0; i < working.size(); i++){
                AlignmentLine curr = working.get(i);
                //Change the number after the greater than sign to change the scaffold overlap threshold
                if((curr.Start() - maxEnd) > 0){
                    if(minStart > 0){
                        if(chr.equals("chr13")){
                            for(AlignmentLine a : aStorage){
                                if(a.Name().equals("7180008663831")){
                                    //System.err.println("hey");
                                }
                            }
                        }
                        this.mergeStore.add(new AlignBed(chr, minStart, maxEnd, aStorage));
                        // reset
                        aStorage = new ArrayList<>();
                    }
                    minStart = curr.Start();
                    maxEnd = curr.End();
                    aStorage.add(curr);
                }else{
                    merger++;
                    if(curr.End() > maxEnd){
                        maxEnd = curr.End();
                    }
                    aStorage.add(curr);
                }
                prev = curr;
            }
            if(minStart >= 0){
                
                this.mergeStore.add(new AlignBed(chr, minStart, maxEnd, aStorage));
            }
            
            

        }
        System.err.println("[ScaffBuilder] Merged " + merger + " Alignment lines!");
    }
    
    // Different logic model that extends overlapping scaffolds to create final scaffold list
    public void seedExtendScaffolds(){
        // Create the scaffold name tree
        int mergecounter = 0, missassembles = 0;
        alignLineDeconvolution();
        
        // Create the alignbed tree for easy reference
        alignTreeGeneration();
                
        // major loop to go through every merged list
        Set<String> used = new HashSet<>();
        for(String scafName : this.sArray.keySet()){ 
            if(scafName.equals("7180008666757")){
                //System.err.println("hey");
            }
            if(used.contains(scafName)){continue;} // Skip if we've already used this one in an other merger
        // 1. Take care of single aligned scaffolds first
            if(alignTree.containsKey(scafName) && alignTree.get(scafName).size() <= 1 && 
                    !(this.misassemblies.contains(scafName)) &&
                    !(mergedSample.contains(scafName))){
                
                SuperScaffold s = new SuperScaffold();
                s.addAlignments(alignTree.get(scafName).get(0).alignments.get(0), scafName);
                alignTree.get(scafName).get(0).SetToUsed();
                this.scaffolds.add(s);
                continue;
            }else{
        // 2. Now, seed extend logic
                // Check if there are overlapping entries
                Set<String> mergeScaffs = null;
                for(Set<String> s : alignStore){
                    if(s.contains(scafName)){
                        mergeScaffs = s;
                        break;
                    }
                }
                ArrayList<AlignmentLine> storage = new ArrayList<>();
                HashMap<String, ArrayList<AlignmentLine>> uniqueList = new HashMap<>();
                
                // Adding all scaffolds to a storage array for later sorting
                // Select the Seed or pick the first element if there isn't a seed
                ArrayList<SuperScaffold> scafStore = new ArrayList<>(); 
                SuperScaffold scaffold = null;
                if(mergeScaffs != null){
                    used.addAll(mergeScaffs);
                    for(String n : mergeScaffs){
                        if(alignTree.containsKey(n) && scaffold == null){
                            scaffold = new SuperScaffold();                            
                            scaffold.addAlignments(alignTree.get(n).get(0).alignments.get(0), alignTree.get(n).get(0).alignments.get(0).Name());
                            // Add all elements from this alignBed
                            for(int x = 1; x < alignTree.get(n).get(0).alignments.size(); x++){
                                //storage.add(alignTree.get(n).get(0).alignments.get(x));
                                if(!uniqueList.containsKey(alignTree.get(n).get(0).alignments.get(x).Chr()))
                                    uniqueList.put(alignTree.get(n).get(0).alignments.get(x).Chr(), new ArrayList<AlignmentLine>());
                                uniqueList.get(alignTree.get(n).get(0).alignments.get(x).Chr()).add(alignTree.get(n).get(0).alignments.get(x));
                            }
                            // Add all elements from any other connected alignBed
                            for(int x = 1; x < alignTree.get(n).size(); x++){
                                for(AlignmentLine l : alignTree.get(n).get(x).alignments){
                                    //storage.add(l);
                                    if(!uniqueList.containsKey(l.Chr()))
                                        uniqueList.put(l.Chr(), new ArrayList<AlignmentLine>());
                                    uniqueList.get(l.Chr()).add(l);
                                }
                            }
                        }else{
                            for(AlignmentLine l : this.sArray.get(n)){
                                //storage.add(l);
                                if(!uniqueList.containsKey(l.Chr()))
                                    uniqueList.put(l.Chr(), new ArrayList<AlignmentLine>());
                                uniqueList.get(l.Chr()).add(l);
                            }
                        }
                    }
                }else{
                    used.add(scafName);
                    if(!this.sArray.containsKey(scafName))
                        System.err.println("Could not find: " + scafName + " in alignmentline factory!");
                    for(AlignmentLine l : this.sArray.get(scafName)){
                        if(scaffold == null){
                            scaffold = new SuperScaffold();
                            scaffold.addAlignments(l, scafName);
                        }else{
                            //storage.add(l);
                            if(!uniqueList.containsKey(l.Chr()))
                                uniqueList.put(l.Chr(), new ArrayList<AlignmentLine>());
                            uniqueList.get(l.Chr()).add(l);
                        }
                    }
                }
                
                ArrayList<String> sortedChrs = utils.SortByChr.ascendingChr(uniqueList.keySet());
                
                for(String chr : sortedChrs){
                    ArrayList<AlignmentLine> aList = uniqueList.get(chr);
                    Collections.sort(aList, new Comparator<AlignmentLine>(){

                        @Override
                        public int compare(AlignmentLine t, AlignmentLine t1) {
                            return t.Start() - t1.Start();
                        }
                        
                    });
                    AlignmentLine last = null;
                    for(AlignmentLine l : aList){
                        
                        if(last == null){
                            if(!scaffold.containsAlign(l))
                                storage.add(l);
                            last = l;
                            continue;
                        }else{
                            if(last.equals(l))
                                continue;
                            if(!scaffold.containsAlign(l))
                                storage.add(l);
                            last = l;
                        }
                    }
                }
                
                scafStore.add(scaffold);
                if(storage.size() == 0){
                    this.scaffolds.addAll(scafStore);
                    continue;
                }
                
                // If genunine difference, or overlapping segment with a missassembly on the two requisite sides, add to the current scaffold list
                // OK, now we do the loop through the alignlist and we add things to the superscaffolds
                // The alignTree holder serves as an excellent extender, but there will be cases where the scaffold cannot be extended
                int storeUsed[] = new int[storage.size()];
                Arrays.fill(storeUsed, 0);
                int num = 0; // counter for the number of used scaffolds

                while(true){                    
                    if(num == storage.size())
                        break; // We've reached the end of the alignment line list. Time to end the loop
                    
                    boolean found = false;
                    // First, try to add to existing seed scaffold
                    for(int x = 0; x < storage.size(); x++){
                        if(storeUsed[x] == 1) continue;
                        
                        AlignmentLine working = storage.get(x);
                        for(SuperScaffold s : scafStore){
                            // Start first
                            if(s.overlaps(working, ovlpBaseNum, genuines)){
                                storeUsed[x] = 1;
                                num++;
                                found = true;
                                break;
                            }
                        }
                        if(found)
                            break; // We extended the seed. Time to start over again in the next loop!
                    }
                    if(found)
                        continue;
                    
                    // If the scaffold cannot be extended, then take the next available alignment and create a new superscaffold
                    for(int x = 0; x < storage.size(); x++){
                        if(storeUsed[x] == 1) continue;
                        SuperScaffold temp = new SuperScaffold();
                        temp.addAlignments(storage.get(x), storage.get(x).Name());
                        if(storage.get(x).Name().equals("7180008666757")){
                            System.err.println("hey");
                        }
                        scafStore.add(temp);
                        storeUsed[x] = 1;
                        num++;
                        break;
                    }
                    
                }        
                // Add scaffolds to class attribute scaffolds
                this.scaffolds.addAll(scafStore);
            }
        }
    }
    
    @Deprecated
    public void chainScaffolds(){
        // 1. Take care of single aligned scaffolds and genuine differences
        // Create the scaffold name tree
        HashMap<String, ArrayList<AlignBed>> alignTree = new HashMap<>();
        
        for(AlignBed a : this.mergeStore){
            for(String name : a.getNames()){
                
                if(!alignTree.containsKey(name)){
                    alignTree.put(name, new ArrayList<AlignBed>());
                }
                alignTree.get(name).add(a);
            }
        }
        int mergecounter = 0, missassembles = 0;
        
        
        
        
        // Go through by scaffold name and count number of alignments associated
        for(String scafName : alignTree.keySet()){
            
            // Scaffolds that only align in one spot and have not been used to create other super scaffolds are binned automatically
            if(alignTree.get(scafName).size() <= 1 
                    && !alignTree.get(scafName).get(0).used 
                    && !this.misassemblies.contains(scafName)){
                SuperScaffold s = new SuperScaffold();
                s.addAlignments(alignTree.get(scafName).get(0).alignments.get(0), scafName);
                alignTree.get(scafName).get(0).SetToUsed();
                this.scaffolds.add(s);
                continue;
            }else if (this.genuines.contains(scafName)){
                // Scaffolds that have multiple mappings are checked to see if they are a genuine difference. If so, all of their alignments are turned into a super scaffold
                boolean found = false;
                for(AlignBed a : this.mergeStore){
                    int total = 0, same = 0, gen = 0;
                    for(AlignmentLine b : a.alignments){
                        total++;
                        if(b.Name().equals(scafName)){
                            same++; gen++;
                        }else if(!b.Name().equals(scafName)){
                            gen++;
                        }
                    }
                    if(same > 1){
                        if(same != total && gen != total){
                            // We have a genuine difference that is together with non-genuine differences; we need to consider them carefully
                            found = true;
                        }
                    }
                }
                if(!found){
                    // Scaffold is a genuine difference, and all alignments should be placed into a superscaffold
                    SuperScaffold s = new SuperScaffold();
                    for(AlignBed a : alignTree.get(scafName)){
                        s.addAlignments(a.alignments, a.getNames());
                        a.SetToUsed();
                    }
                    this.scaffolds.add(s);
                }
            }
        }
        
        // 2. Now take care of the remaining batch that have misassemblies. 
        // I need to scroll through by the alignment coordinates, and then conditionally merge
        // I need an array containing just the misassemblies/oddities for further processing
        
        for(String scafName : alignTree.keySet()){
            if(scafName.equals("7180008670784")){
                System.out.println("hey");
            }

            // Collect the misassemblies that I need to work on
            // So, I need to check merger across scaffold names and chromosome coordinates here
            ArrayList<AlignmentLine> misassembleStore = new ArrayList<>();
            for(AlignBed a : alignTree.get(scafName)){
                if(a.used) continue;
                else{
                    Set<String> names = new HashSet<>();
                    for(AlignmentLine v : a.alignments){
                        if(!v.Name().equals(scafName) && !v.wasUsed){
                            names.add(v.Name());
                        }
                    }
                    if(!names.isEmpty()){
                        for(String n : names){
                            for(AlignBed b : alignTree.get(n)){
                                if(!b.used){
                                    for(AlignmentLine c : b.alignments){
                                        if(!c.wasUsed){
                                            misassembleStore.add(c);
                                        }
                                    }
                                    b.SetToUsed();
                                }
                            }
                        }
                    }
                }
                for(AlignmentLine b : a.alignments){
                    if(!b.wasUsed){
                        misassembleStore.add(b);
                    }
                }
               a.SetToUsed();
            }
            
            
            Collections.sort(misassembleStore, new Comparator<AlignmentLine>(){

                @Override
                public int compare(AlignmentLine t, AlignmentLine t1) {
                    if(t.Name().equals(t1.Name())){
                        int myStart = (t.reversedCoords())? t.ContigEnd() : t.ContigStart();
                        int yourStart = (t1.reversedCoords())? t1.ContigEnd() : t1.ContigStart();
                        return yourStart - myStart;
                    }else{
                        if(!t.Chr().equals(t1.Chr())){
                            int mychr = utils.SortByChr.retVal(t.Chr());
                            int yourchr = utils.SortByChr.retVal(t1.Chr());
                            return mychr - yourchr;
                        }
                        return t.Start() - t1.Start();
                    }
                }
                
            });
            ArrayList<SuperScaffold> scaffoldstore = new ArrayList<>();
            // Go through each alignmentLine in the scaffold
            for(int x = 0; x < misassembleStore.size(); x++){
                // scaffold "start" or "end" misassemblies are based on scaffold alignment coordinates
                // end misassemblies are based on the "contigEnd" attribute; vice versa for start

                AlignmentLine current = misassembleStore.get(x);


                if(scaffoldstore.isEmpty()){
                    // Add scaffold to the array if there are  none in there already
                    SuperScaffold s = new SuperScaffold();
                    s.addAlignments(current, current.Name());
                    scaffoldstore.add(s);
                }else{
                    // There are other scaffolds, check to see if we need to update existing ones or create a new one
                    boolean found = false;
                    for(SuperScaffold s : scaffoldstore){
                        for(AlignmentLine b : s.getOutScaffoldListSegments()){
                            //Check to see if it is a genuine difference and if it belongs to the same scaffold number
                            if(current.Name().equals(b.Name())){
                                s.addAlignments(current, current.Name());
                                found = true;
                                mergecounter++;
                                break;
                            }else{
                                // Now we get into the complicated logic
                                int curSCstart = (current.reversedCoords())? current.ContigEnd(): current.ContigStart();
                                boolean curSMstart = (current.reversedCoords())? current.MissEnd() : current.MissStart();
                                boolean curSMend = (current.reversedCoords())? current.MissStart() : current.MissEnd();
                                int bSCstart = (b.reversedCoords())? b.ContigEnd(): b.ContigStart();
                                boolean bSMstart = (b.reversedCoords())? b.MissEnd() : b.MissStart();
                                boolean bSMend = (b.reversedCoords())? b.MissStart() : b.MissEnd();

                                boolean sameName = current.Name().equals(b.Name());
                                boolean sameChr = current.Chr().equals(b.Chr());


                                if((current.Start() < b.Start() && current.MissEnd() && b.MissStart() && !sameName && sameChr) ||
                                        (b.Start() < current.Start() && current.MissStart() && b.MissEnd() && !sameName && sameChr)){
                                    // There was a break between the scaffolds and they don't share the same name
                                    // We join them because we believe that the misassembly indicates that they should be stuck together
                                    s.addAlignments(current, current.Name());
                                    found = true;
                                    mergecounter++;
                                    break;
                                }else if((sameName) &&
                                        ((curSCstart < bSCstart && !curSMend && !bSMstart) ||
                                        (bSCstart < curSCstart && !curSMstart && !bSMend))){
                                    // The scaffolds are named the same, and there are no misassemblies between them
                                    s.addAlignments(current, current.Name());
                                    found = true;
                                    mergecounter++;
                                    break;
                                } 
                            }
                        }

                        if(found) break; // Just to avoid adding the same alignmentline to multiple scaffolds
                    }
                    if(!found){
                        // The scaffold didn't meet our criteria, so we create a new superscaffold for it and leave it be
                        SuperScaffold s = new SuperScaffold();
                        s.addAlignments(current, current.Name());
                        scaffoldstore.add(s);
                        missassembles++;
                    }
                }
            }
            this.scaffolds.addAll(scaffoldstore);
            
        }
        
        System.err.println("[ScaffBuilder] Identified " + mergecounter + " mergers and " + missassembles + " misassembled partitions!");
        
    }
    
    
    // Second step! Break apart the misassemblies and rejoin if they overlap with another segment
    @Deprecated
    public void linkUpScaffolds(){
        // Create the scaffold name tree
        HashMap<String, ArrayList<AlignBed>> alignTree = new HashMap<>();
        
        for(AlignBed a : this.mergeStore){
            for(String name : a.getNames()){
                if(!alignTree.containsKey(name)){
                    alignTree.put(name, new ArrayList<AlignBed>());
                }
                alignTree.get(name).add(a);
            }
        }
        int mergecounter = 0, missassembles = 0;
        
        // Go through by scaffold name and count number of alignments associated
        for(String scafName : alignTree.keySet()){
            // Scaffolds that only align in one spot and have not been used to create other super scaffolds are binned automatically
            if(alignTree.get(scafName).size() <= 1 
                    && !alignTree.get(scafName).get(0).used 
                    && !this.misassemblies.contains(scafName)){
                SuperScaffold s = new SuperScaffold();
                s.addAlignments(alignTree.get(scafName).get(0).alignments.get(0), scafName);
                alignTree.get(scafName).get(0).used = true;
                this.scaffolds.add(s);
                continue;
            }else{
                boolean used = false;
                for(AlignBed a : alignTree.get(scafName)){
                    if(a.used){
                        used = true;
                        break;
                    }
                }
                if(used){
                    // This scaffold has been incorporated into another superscaffold
                    mergecounter++;
                    continue;
                }
                // Scaffolds that have multiple mappings are checked to see if they are a genuine difference. If so, all of their alignments are turned into a super scaffold
                if(this.genuines.contains(scafName)){
                    // Scaffold is a genuine difference, and all alignments should be placed into a superscaffold
                    SuperScaffold s = new SuperScaffold();
                    for(AlignBed a : alignTree.get(scafName)){
                        s.addAlignments(a.alignments, a.getNames());
                        a.used = true;
                    }
                    this.scaffolds.add(s);
                }else if(this.misassemblies.contains(scafName)){
                    // Scaffold is a missassembly and needs to be chopped up
                    ArrayList<AlignmentLine> aligns = new ArrayList<>();
                    ArrayList<String> names = new ArrayList<>();
                    if(scafName.equals("7180008666595")){
                        System.out.println("hey!");
                    }
                    for(AlignBed array : alignTree.get(scafName)){
                        array.used = true;
                        Collections.sort(array.alignments);
                        for(AlignmentLine l : array.alignments){
                            if(l.MissStart()){
                                // Missassembly at the start of this scaffold; break and make a new one
                                if(!aligns.isEmpty()){
                                    SuperScaffold s = new SuperScaffold();
                                    s.addAlignments(aligns, names);
                                    this.scaffolds.add(s);
                                    aligns = new ArrayList<>();
                                    names = new ArrayList<>();
                                    missassembles++;
                                }
                                
                                aligns.add(l);
                                names.add(l.Name());
                            }else if(l.MissEnd()){
                                // Missassembly at the end of this scaffold
                                aligns.add(l);
                                names.add(l.Name());
                                
                                SuperScaffold s = new SuperScaffold();
                                s.addAlignments(aligns, names);
                                this.scaffolds.add(s);
                                aligns = new ArrayList<>();
                                names = new ArrayList<>();
                                missassembles++;
                            }else{
                                aligns.add(l);
                                names.add(l.Name());
                            }
                        }
                        
                    }
                    if(!aligns.isEmpty()){
                        SuperScaffold s = new SuperScaffold();
                        s.addAlignments(aligns, names);
                        this.scaffolds.add(s);
                        missassembles++;
                    }
                }else{
                    // Wait until the next iteration to create superscaffolds
                }
            }
        }
        
        // Second iteration to catch all remaining scaffolds that have not been used yet
        for(String scafName : alignTree.keySet()){
            boolean used = false;
            for(AlignBed a : alignTree.get(scafName)){
                if(a.used){
                    used = true;
                    break;
                }
            }
            if(used){
                // This scaffold has already been incorporated into another superscaffold
                continue;
            }else{
                // This is a remaining scaffold and it needs to be added to the pile
                SuperScaffold s = new SuperScaffold();
                for(AlignBed a : alignTree.get(scafName)){
                    s.addAlignments(a.alignments, a.getNames());
                    a.used = true;
                }
                this.scaffolds.add(s);
            }
        }
        System.err.println("[ScaffBuilder] Identified " + mergecounter + " mergers and " + missassembles + " misassembled partitions!");
        
    }
    
    public int ovCount (int start1, int end1, int start2, int end2){
        return soonest(end1, end2) - latest(start1, start2);
    }
    private int soonest (int a, int b){
        return (a >= b)? b : a;
    }
    private int latest (int a, int b){
        return (a >= b)? a : b;
    }

    private void alignLineDeconvolution() {
        // This will just be brute force association. I don't know a more elegant way
        for(AlignBed a : this.mergeStore){
            if(alignStore.isEmpty()){
                Set<String> temp = new HashSet<>();
                temp.addAll(a.getNames());
                alignStore.add(temp);
                mergedSample.addAll(a.getNames());
                continue;
            }
            mergedSample.addAll(a.getNames());
            boolean found = false;
            for(Set<String> aligns : alignStore){                
                for(String name : a.getNames()){
                    if(aligns.contains(name)){
                        aligns.addAll(a.getNames());
                        found = true;
                    }
                }
            }
            
            if(!found){
                Set<String> temp = new HashSet<>();
                temp.addAll(a.getNames());
                alignStore.add(temp);
            }
        }
    }
    private void alignLineHasher(AlignmentLineFactory fact){
        for(String chr : fact.getAllAligns().keySet()){
            for(AlignmentLine l : fact.getAlignFromChr(chr)){
                if(!this.sArray.containsKey(l.Name()))
                    this.sArray.put(l.Name(), new ArrayList<AlignmentLine>());
                this.sArray.get(l.Name()).add(l);
            }
        }
    }
    private void alignTreeGeneration() {
        for(AlignBed a : this.mergeStore){
            for(String name : a.getNames()){
                
                if(!alignTree.containsKey(name)){
                    alignTree.put(name, new ArrayList<AlignBed>());
                }
                alignTree.get(name).add(a);
            }
        }
    }
    
    private class AlignBed extends BedAbstract{
        public ArrayList<AlignmentLine> alignments = new ArrayList<>();
        public boolean used = false;
        
        public AlignBed(String chr, int start, int end, ArrayList<AlignmentLine> aligns){
            this.chr = chr;
            this.start = start;
            this.end = end;
            this.alignments = aligns;
        }
        
        public AlignBed(ArrayList<AlignBed> aligns){
            Collections.sort(aligns);
            for(AlignBed v : aligns){
                this.alignments.addAll(v.alignments);
            }
            this.chr = aligns.get(0).chr;
            this.start = aligns.get(0).start;
            this.end = aligns.get(0).end;
        }
        
        public AlignBed(String chr, int start, int end, AlignmentLine a){
            this.chr = chr;
            this.start = start;
            this.end = end;
            this.alignments.add(a);
        }
        
        @Override
        public int compareTo(BedAbstract t) {
            if(!this.Chr().equals(t.Chr())){
                int mychr = utils.SortByChr.retVal(chr);
                int yourchr = utils.SortByChr.retVal(t.Chr());
                return yourchr - mychr;
            }
            return t.Start() - this.Start();
        }
        
        public ArrayList<String> getNames(){
            ArrayList<String> names = new ArrayList<>();
            for(AlignmentLine l : this.alignments){
                names.add(l.Name());
            }
            return names;
        }
        
        public void SetToUsed(){
            for(AlignmentLine l : this.alignments){
                l.wasUsed = true;
            }
            this.used = true;
        }
    }
    
    public void printOutScaffoldPlan(String outputName){
        try(BufferedWriter output = Files.newBufferedWriter(Paths.get(outputName), Charset.defaultCharset())){
            Collections.sort(scaffolds);
            ArrayList<Integer> scafSizes = new ArrayList<>();
            for(int x = 0; x < scaffolds.size(); x++){
                SuperScaffold s = scaffolds.get(x);
                scafSizes.add(s.length);
                ArrayList<AlignmentLine> list = s.getOutScaffoldListSegments();
                for(AlignmentLine l : list){
                    // Output format:
                    // SuperScaffold SubScaffoldName SScaffoldStart SScaffoldEnd bChr bStart bEnd
                    output.write("Scaffold" + (x+1) + "\t" + l.Name() + "\t" + l.ContigStart() + "\t"
                            + l.ContigEnd() + "\t" + l.Chr() + "\t" + l.Start() + "\t" + l.End() + "\t" + 
                            l.MissStart() + "\t" + l.MissEnd() + "\t" + System.lineSeparator());
                }
            }
            int n50 = IMedian(scafSizes);
            System.out.println("Scaffold N50 size = " + String.format("%.2f kb", (float) n50 / 1000));
        }catch(IOException ex){
            ex.printStackTrace();
        }
    }
    
    public int IMedian(List<Integer> list){       
        if(list.size() == 1){
            return list.get(0);
        }else if(list.isEmpty()){
            return 0;
        }else{
            Collections.sort(list);
            if(list.size() % 2 == 0){
                return list.get(list.size() / 2);
            }else{
                return list.get( (int)(list.size() / 2) + 1);
            }
        }        
    }
}

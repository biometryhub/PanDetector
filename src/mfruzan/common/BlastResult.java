/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package mfruzan.common;

import htsjdk.samtools.reference.FastaSequenceFile;
import htsjdk.samtools.reference.FastaSequenceIndex;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequence;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.Reader;
import java.io.Writer;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.ListIterator;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.zip.GZIPInputStream;
import mfruzan.algorithm.BlastHSP;

/**
 *
 * @author mario.fruzangohar
 * This class process Blast result file in table format(6 and 7) and also lastZ file MAF
 */
public class BlastResult {
    BufferedReader reader = null;
    private String fileName = null;
    byte format = 0; 
    public BlastResult(String fileName, byte format){
        this.fileName = fileName;
        this.format = format;
        try{
            if (fileName.endsWith(".gz")){
                InputStream filein = new FileInputStream(fileName);
                InputStream gzipin = new GZIPInputStream(filein);
                Reader ir = new InputStreamReader(gzipin);
                this.reader = new BufferedReader(ir);
            }else
                this.reader = new BufferedReader(new FileReader(this.fileName));
        } catch (Exception ex) {
            System.out.println(ex.getMessage());
        }
    }
    public BlastResult(){
        
    }

    
    
    private void insertIntoBlastHSPList2(LinkedList<BlastHSP> list, BlastHSP aHSP){
        
        aHSP.normalizeOrientation(); // if orientation is +/- then swap sstart and send
        
        boolean added = false;
        for(int i=0;i<list.size();i++){
            BlastHSP hsp = list.get(i);
            if(aHSP.s_start < hsp.s_start){
                list.add(i, aHSP);
                added = true;
                break;
            }
        }
        
        if (!added)// then it means that aHSP comese after all the elements of the list
            list.add(aHSP);
    }
    
    public static boolean insertIntoBlastHSPList(LinkedList<BlastHSP> list, BlastHSP aHSP){
        int i=0;
        for(; i<list.size(); i++){
            BlastHSP hsp = list.get(i);
            if (i==0){
                if (aHSP.comesBefore(hsp, 12)){
                    list.add(0, aHSP);
                    return true;
                }
            }else{
                // i>0, then aHSP should be between element i-1 and i to be able to be inserted into i
                if (aHSP.comesBefore(hsp, 12) && aHSP.comesAfter(list.get(i-1), 12)){
                    list.add(i, aHSP);
                    return true;
                }
            }
        }
        // we get here for 3 reasons
        //1: list is empty
        if (list.size()==0){
            list.add(aHSP);
            return true;
        }
        //2: aHSP comes after last element of the list
        if (aHSP.comesAfter(list.get(list.size()-1), 12)){
            list.add(aHSP);// add it 
            return true;
        }
        //3: aHSP has significant overlap with existing HSPs in the list
        
        return false;
    }

  
}

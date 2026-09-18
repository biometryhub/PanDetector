/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package mfruzan.algorithm;

import mfruzan.common.Helper;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.NavigableSet;
import java.util.Set;
import java.util.TreeMap;

/**
 *
 * @author mario.fruzangohar
 * This class encapsulate blast high-scoring segment pairs (HSP), it also can encapsulate binary MAF record 
 */
public class BlastHSP {
    public String q_contig = null;
    public int q_start = 0;  //q_start is always less than q_end, but not for the subject
    public int q_end = 0;
    public int q_len = 0;
    
    public String s_contig = null;
    public int s_start = 0;
    public int s_end = 0;
    public int s_len = 0;
    
    public int alignment_length = 0;  // longer the better
    public double identity_percent = 0; // a number between 0-1
    public double e_value = 1;
    public double score = 0;
    
    // Alignment block can contain - that denotes gaps in the alignment. MAF file format provide us with this type of information.
    public String q_alignment_block = null;
    public String s_alignment_block = null;
    public boolean q_pos_strand = true; // used in maf format
    public boolean s_pos_strand = true; // used in maf format
    public int q_contig_size = 0; // used in maf format
    public int s_contig_size = 0; // used in maf format
    //calculated field (is used in MAF file format)
    public byte orientation = 1; // 1:same (+/+ or -/-) 0:reverse (+/- or -/+); this is to see how query and subject are aligned in relation to each other
    
    public String scoreLine = null;
    public String subjectLine = null;
    public String queryLine = null;
    // If MAF file is multi-line to represent Multiple sequence alignment then this variable is used
    public List<MAFSLine> mafLines = null;
    
    
    // constructor
    public BlastHSP(){
       
    }
    @Override
    public BlastHSP clone(){
        BlastHSP ahsp = new BlastHSP();
        ahsp.q_contig = this.q_contig;
        ahsp.q_start = this.q_start;  //q_start is always less than q_end, but not for the subject
        ahsp.q_end = this.q_end;
        ahsp.q_len = this.q_len;

        ahsp.s_contig = this.s_contig;
        ahsp.s_start = this.s_start;
        ahsp.s_end = this.s_end;
        ahsp.s_len = this.s_len;

        ahsp.alignment_length = this.alignment_length;  // longer the better
        ahsp.identity_percent = this.identity_percent; // a number between 0-1
        ahsp.e_value = this.e_value;
        ahsp.score = this.score;

        ahsp.q_alignment_block = this.q_alignment_block;
        ahsp.s_alignment_block = this.s_alignment_block;
        ahsp.q_pos_strand = this.q_pos_strand; 
        ahsp.s_pos_strand = this.s_pos_strand; // used in maf format
        ahsp.q_contig_size = this.q_contig_size; // used in maf format
        ahsp.s_contig_size = this.s_contig_size; // used in maf format
        ahsp.orientation = this.orientation; // 1:same (+/+ or -/-) 0:reverse (+/- or -/+); this is to see how query and subject are aligned in relation to each other

        ahsp.scoreLine = this.scoreLine;
        ahsp.subjectLine = this.subjectLine;
        ahsp.queryLine = this.queryLine;
        ahsp.mafLines = this.mafLines;
        return ahsp;
    }
    public void print4debug(){
       
        System.out.println(String.format("[q_contig:%s,q_start:%d,q_end:%d,q_pos_strand:%b][s_contig:%s,s_start:%d,s_end:%d,s_pos_strand:%b]", q_contig, q_start, q_end, q_pos_strand, s_contig, s_start, s_end, s_pos_strand));
    }
    
    public void printQueryAsMAF(){
         //convert all the cordinates to MAF standard
           int q_zero_start = 0;
           if(this.q_pos_strand){
               q_zero_start = this.q_start - 1;
           }else{
               q_zero_start = this.q_contig_size - this.q_end;
           }
           // now write s line
          System.out.println(String.format("%s\t%d\t%d\t%s\t%d", this.q_contig, q_zero_start, this.q_len, this.q_pos_strand?"+":"-", this.q_contig_size));

    }
    public void printSubjectAsMAF(){
         //convert all the cordinates to MAF standard
           int s_zero_start = 0;
           if(this.s_pos_strand){
               s_zero_start = this.s_start - 1;
           }else{
               s_zero_start = this.s_contig_size - this.s_end;
           }
           // now write s line
          System.out.println(String.format("%s\t%d\t%d\t%s\t%d", this.s_contig, s_zero_start, this.s_len, this.s_pos_strand?"+":"-", this.s_contig_size));

    }
    
    public BlastHSP(String qcontig, int qs, int qe, String scontig, int ss, int se, int len, double percent, double eval){
        q_contig = qcontig;
        q_start = qs;
        q_end = qe;
        s_contig = scontig;
        s_start = ss;
        s_end = se;
        alignment_length = len;
        identity_percent = percent;
        e_value = eval;
        if (s_start>s_end)
            orientation = -1;
    }
    
    public String getSubjectPart(int start, int end, String query_pattern) throws Exception{
        //returns part of subject corresponding query start and end (both 1-based)
        int start_idx = 0;
        int end_idx = 0;
        int pointer = 0;
        for(int i=0; i<q_alignment_block.length(); i++){
            if(q_alignment_block.charAt(i) != '-')
                pointer++; 
            if (pointer==start)
                start_idx = i;
            if (pointer==end){
                end_idx = i;
                break;
            }
              
            
        }
 
        String initial_block = s_alignment_block.substring(start_idx, end_idx+1);
        
        // If it is longer than query_pattern then raise the exception:
        if (initial_block.length() > query_pattern.length()){
            System.out.println("extracted subject is longer than initial query pattern " + initial_block.length() + " > " + query_pattern.length());
            System.out.println("query pattern:" + query_pattern);
            System.out.println("subject part:" + initial_block);
            
            System.out.println(this.toString());
        }
        // see if it is shorter than query_pattern then inject - from query_pattern until length get equal
        if (initial_block.length() < query_pattern.length()){
            StringBuffer buff = new StringBuffer();
            int diff = query_pattern.length() - initial_block.length();
            int j=0; // j is index of initial_block
            for (int i = 0; i<query_pattern.length(); i++){
                if (query_pattern.charAt(i)=='-' && diff>0){
                    buff.append('-');
                    diff--;
                }else{
                    buff.append(initial_block.charAt(j));
                    j++;
                }
            }
            
            if(buff.length() != query_pattern.length())
                throw new Exception("Still extracted subject is not equal length of expected query ");
            
            return buff.toString();
        }//if  
        return initial_block.toString();
    }
    
    public void alignSubjectPart(int start, int end, int index, List<StringBuffer> fragments) throws Exception{
        if (index == 3 &&  s_contig.equals("NQIK01000004.1"))
            index = 3;
       
        int start_idx = 0; // 0-based start of subject equivalent of query
        int end_idx = 0;  // 0-based end of subject equivalent of query
        int pointer = 0; // 0-based pointer at query string (q-qlignment-block minus its -)
        Map<Integer, Integer> idxs = new TreeMap(); // Map of postions and number of times where - observed (TreeMap orderd by Key)
        TreeMap<Integer, Integer> adjustedIdxs = null;
        //StringBuffer buff = new StringBuffer();//this buff wil hold extracted_subject
        for(int i=0; i<q_alignment_block.length(); i++){
            if(q_alignment_block.charAt(i) != '-')
                pointer++; 
            //else
            //    idxs.add(i);
            if (pointer==start)
                start_idx = i;
            if (pointer==end){
                end_idx = i;
                break;
            }        
            //when we get here it means we have not reached to pointer end, but we may be passed start, if we passed start and current character is - then put its relative 0-based index into inxs
            if (pointer>=start && q_alignment_block.charAt(i) == '-')
                new MapUtil().increment(idxs, pointer-start + 1); 
                
        }
 
        String extracted_subject = s_alignment_block.substring(start_idx, end_idx+1);
        if (extracted_subject.length() == fragments.get(0).length() ){     
            fragments.get(index).append(extracted_subject);
            return;
        } 
        
        String  updatedQuery = fragments.get(0).toString();
        if (extracted_subject.length() != fragments.get(0).length() && idxs.size()>0){
            adjustedIdxs = adjustIndex(fragments.get(0).toString(), idxs);
            StringBuffer buff = new StringBuffer(fragments.get(0));
            for (int i=0; i<index; i++){
                for (Entry<Integer, Integer> ent : adjustedIdxs.descendingMap().entrySet()){
                    if (ent.getValue().intValue()>0){
                        for (int j=0; j<ent.getValue().intValue(); j++){                        
                           fragments.get(i).insert(ent.getKey().intValue(), '-');
                        }
                    }else if (ent.getValue().intValue()<0 ){
                        if (i==0)
                            buff.delete(ent.getKey().intValue()+ent.getValue().intValue(), ent.getKey().intValue());
                                                
                    }
                }
            }
            
             updatedQuery = buff.toString();
            
        }
        if (fragments.get(0).length() == extracted_subject.length()){
            fragments.get(index).append(extracted_subject);
            return;
        }
        
        if (extracted_subject.length() < fragments.get(0).length()){
           
            StringBuffer buff = new StringBuffer();
            int j=0; // j is index of extracted_subject
            for (int i = 0; i<updatedQuery.length(); i++){
                if (updatedQuery.charAt(i)=='-' /*&& extracted_subject.charAt(j)!='-'*/){
                    buff.append('-');
                }else{
                    if (j<extracted_subject.length())
                      buff.append(extracted_subject.charAt(j++));
                    
                }
            }
            if (j<extracted_subject.length())
                buff.append(extracted_subject.substring(j));
            fragments.get(index).append(buff);
            if(buff.length() != fragments.get(0).length())
                throw new Exception("Still extracted subject is not equal length to expected query at level " + index + " subject " +s_contig + " Start " + s_start + " End " + s_end );                        
            if (Helper.getSimilarity(fragments.get(0).toString(), buff.toString()) <0.9)
                System.out.println("extension for subject wrong at level " + index + " subject " +s_contig + " Start " + s_start + " End " + s_end);
            
        }
      
        
    }
   
    
    public TreeMap<Integer, Integer> adjustIndex(String seq, Map<Integer, Integer> idxs){
        TreeMap<Integer, Integer> out = new TreeMap();
        
        int pointer = -1; // is original(gapless) index    of seq     
        int gappos = -1;
        int count = -1;
        Entry<Integer, Integer> ent = null;
        Iterator<Entry<Integer,Integer>> it = idxs.entrySet().iterator();
        if (it.hasNext()){
            ent = it.next();
            gappos = ent.getKey();
            count = ent.getValue();
                    
        }
        for (int i=0; i<seq.length();i++){
            if (seq.charAt(i)!='-')
                pointer++;

            if (pointer == gappos){
                int dash_cnt = 0;
                if (i>0)
                    for (int j=i-1; j>=0; j--)
                        if (seq.charAt(j)=='-')
                            dash_cnt++;
                        else
                            break;
                
                int remainder = count - dash_cnt;// remainder can be positive, zero or negative, each has its own treatment (all are needed to be returned for the algorithm to work)                
                 out.put(i, remainder);
                if (it.hasNext()){
                    ent = it.next();
                    gappos = ent.getKey();
                    count = ent.getValue();                    
                }else
                    break;
            }
                
        }
        return out;
    }
    public void normalizeOrientation(){
        //if HSP is +/- then send<sstart, then we just swap those to be normal
        if (this.s_start>this.s_end){
            int temp = this.s_start;
            this.s_start = this.s_end;
            this.s_end = temp;
        }
    }
    public boolean distinctSubject(BlastHSP hsp){
        // to see if hsp's subject is distinct from this HSP
        if (!hsp.s_contig.equals(this.s_contig))
            return true;
        int hsp_left = 0;
        int hsp_right = 0;
        if (hsp.s_start<hsp.s_end){
            hsp_left = hsp.s_start;
            hsp_right = hsp.s_end;
        }else{
            hsp_left = hsp.s_end;
            hsp_right = hsp.s_start;            
        }
        
        int s_left = 0;
        int s_right = 0;
        if (this.s_start<this.s_end){
            s_left = s_start;
            s_right = s_end;
        }else{
            s_left = s_end;
            s_right = s_start;
        }
        //---------------------------------------
        if (hsp_left > s_right  || hsp_right<s_left)
            return true;
        
        return false;
    }
    
    public boolean subjetConnected(BlastHSP hsp , int overlap){
        // returns true if this hsp' subject continues the other one and also 2 hsp have separate query regions
        if (!hsp.s_contig.equals(this.s_contig))
            return false;

        if(hsp.q_contig.equals("k99_16175") && (hsp.q_start==35 || hsp.q_start==1114))
            System.out.println();
        if(!distinctQuery(hsp,0))
            return false;
        //first fix start and end of the subject
        int hsp_left = 0;
        int hsp_right = 0;
        if (hsp.s_start<hsp.s_end){
            hsp_left = hsp.s_start;
            hsp_right = hsp.s_end;
        }else{
            hsp_left = hsp.s_end;
            hsp_right = hsp.s_start;            
        }
        
        int s_left = 0;
        int s_right = 0;
        if (this.s_start<this.s_end){
            s_left = s_start;
            s_right = s_end;
        }else{
            s_left = s_end;
            s_right = s_start;
        }
        if (hsp_left<=s_left &&  Math.abs(hsp_right - s_left) <= overlap)
            return true;
        
        if(hsp_left>=s_left && Math.abs(hsp_left - s_right) <= overlap)
            return true;
        
        return false;
    }
    
    public void trimQuery(BlastHSP aHsp, boolean log){
        
        if (!aHsp.q_contig.equals(this.q_contig))
            return ;
        if (aHsp.q_start > this.q_end  || aHsp.q_end<this.q_start)
            return ;
        if (this.q_start >= aHsp.q_start   && this.q_end <= aHsp.q_end ){
            System.out.println("Info:this covered by aHsp");
            this.q_len = 0;
            return ;
            
        }else if (aHsp.q_start >= this.q_start   && aHsp.q_end <= this.q_end){
            System.out.println("Warning:aHsp covered by this");
            this.q_len = 0;
            return ;
        }
        
        try{
            if( aHsp.q_start<=this.q_end && aHsp.q_start>this.q_start){
                int trim_len = this.q_end - aHsp.q_start + 1;
                this.q_end -= trim_len;// q_end may change in the second loop
                int expected_q_len = this.q_len - trim_len;

    
                int i=0;
                int q_pos = this.q_start; // query ref start pos
                int s_pos = this.s_start; // subject ref start pos
                for (; i<this.q_alignment_block.length(); i++){
                    if (this.q_alignment_block.charAt(i)!='-')
                        q_pos++;
                    if (this.s_alignment_block.charAt(i)!='-')
                        s_pos++;

                    if (q_pos > this.q_end) // remember that = will not work
                        break;

                }
                q_pos--;
                for(;i>0; i--){
                    if (this.q_alignment_block.charAt(i)!='-' && this.s_alignment_block.charAt(i)!='-')
                        break;
                    if (this.q_alignment_block.charAt(i)!='-')
                        q_pos--;

                    if (this.s_alignment_block.charAt(i)!='-')
                        s_pos--;
                }
                this.q_end = q_pos;
                this.s_end = s_pos;
                this.alignment_length = i+1;
                this.q_len = this.q_end - this.q_start + 1;
                this.s_len = this.s_end - this.s_start + 1;
                this.q_alignment_block = this.q_alignment_block.substring(0, i+1);
                this.s_alignment_block = this.s_alignment_block.substring(0, i+1);
            }else if( aHsp.q_end>=this.q_start && aHsp.q_end<this.q_end){
                if (log){
                    System.out.println(this.s_alignment_block);
                    System.out.println(this.q_alignment_block);
                }
                int old_q_start = this.q_start;
                int trim_len =  aHsp.q_end  -this.q_start  + 1;
                this.q_start += trim_len;
                int expected_q_len = this.q_len - trim_len;
                int i=0;
                int q_pos = old_q_start; // query ref start pos
                int s_pos = this.s_start; // subject ref start pos
                for (; i<this.q_alignment_block.length(); i++){
                    if (this.q_alignment_block.charAt(i)!='-')
                        q_pos++;
                    if (this.s_alignment_block.charAt(i)!='-')
                        s_pos++;

                    if (q_pos > this.q_start) // remember that = will not work
                        break;

                }
                if (this.s_alignment_block.charAt(i)!='-')
                 s_pos--;
                for(;i<this.q_alignment_block.length(); i++){
                    if (this.q_alignment_block.charAt(i)!='-' && this.s_alignment_block.charAt(i)!='-')
                        break;

                    if (this.q_alignment_block.charAt(i)!='-')
                        q_pos++;

                    if (this.q_alignment_block.charAt(i)!='-')
                        s_pos++;                
                }
                this.q_start = q_pos;
                this.s_start = s_pos;
                this.alignment_length = this.alignment_length - i;
                this.q_len = this.q_end - this.q_start + 1;
                
                this.s_len = this.s_end - this.s_start + 1;
                this.q_alignment_block = this.q_alignment_block.substring(i);
                this.s_alignment_block = this.s_alignment_block.substring(i);
                if (log){
                    System.out.println(this.s_alignment_block);
                    System.out.println(this.q_alignment_block);
                }
                
            }
        
        }catch(Exception ex){
            System.out.println(ex.getMessage());
        }
        
    }
    public void trimQuery2(BlastHSP aHsp, boolean log){
        
        if (!aHsp.q_contig.equals(this.q_contig))
            return ;
        if (aHsp.q_start > this.q_end  || aHsp.q_end<this.q_start)
            return ;
        if (this.q_start >= aHsp.q_start   && this.q_end <= aHsp.q_end ){
            System.out.println("Info:this covered by aHsp");
            this.q_len = 0;
            return ;
            
        }else if (aHsp.q_start >= this.q_start   && aHsp.q_end <= this.q_end){
            // see if there is opposite situation aHsp covered by this
            System.out.println("Warning:aHsp covered by this");
            this.q_len = 0;
            return ;
        }
        
        try{
            if( aHsp.q_start<=this.q_end && aHsp.q_start>this.q_start){
                int trim_len = this.q_end - aHsp.q_start + 1;
            if (!this.q_pos_strand){
                this.q_alignment_block = Helper.getRC(this.q_alignment_block);
                this.s_alignment_block = Helper.getRC(this.s_alignment_block);
            }

    
                int i=this.q_alignment_block.length()-1;
                int nongap_q_bases = 0;
                int nongap_s_bases = 0;
                for (; i>=0; i--){
                    if (nongap_q_bases >= trim_len && this.q_alignment_block.charAt(i)!='-' && this.s_alignment_block.charAt(i)!='-') 
                        break;
                    
                    if (this.q_alignment_block.charAt(i)!='-')
                        nongap_q_bases++;
                    if (this.s_alignment_block.charAt(i)!='-')
                        nongap_s_bases++;

                }
                 this.alignment_length = i+1;
                this.q_alignment_block = this.q_alignment_block.substring(0,i+1);
                this.s_alignment_block = this.s_alignment_block.substring(0,i+1);
                this.alignment_length = this.alignment_length - i;
                if (!this.q_pos_strand){
                    // brought it back to original orientation 
                    this.q_alignment_block = Helper.getRC(this.q_alignment_block);
                    this.s_alignment_block = Helper.getRC(this.s_alignment_block);
                }
                this.q_end -= nongap_q_bases;               
                this.q_len = this.q_end - this.q_start + 1;                
                if (this.q_pos_strand){
                    this.s_end -= nongap_s_bases; 
                }else {
                    this.s_start += nongap_s_bases; 
                }
                this.s_len = this.s_end - this.s_start + 1;
            }else if( aHsp.q_end>=this.q_start && aHsp.q_end<this.q_end){
           
               
                int trim_len =  aHsp.q_end  -this.q_start  + 1;
                if (!this.q_pos_strand){
                    this.q_alignment_block = Helper.getRC(this.q_alignment_block);
                    this.s_alignment_block = Helper.getRC(this.s_alignment_block);
                }
                int i=0;

 
                int nongap_q_bases = 0;
                int nongap_s_bases = 0;
                for (; i<this.q_alignment_block.length(); i++){

                    if (nongap_q_bases >= trim_len && this.q_alignment_block.charAt(i)!='-' && this.s_alignment_block.charAt(i)!='-') 
                        break;

                    if (this.q_alignment_block.charAt(i)!='-')
                        nongap_q_bases++;
                    if (this.s_alignment_block.charAt(i)!='-')
                        nongap_s_bases++;
                }
                this.q_alignment_block = this.q_alignment_block.substring(i);
                this.s_alignment_block = this.s_alignment_block.substring(i);
                 this.alignment_length = this.alignment_length - i;
                if (!this.q_pos_strand){
                    this.q_alignment_block = Helper.getRC(this.q_alignment_block);
                    this.s_alignment_block = Helper.getRC(this.s_alignment_block);
                }
                
                this.q_start += nongap_q_bases;               
                this.q_len = this.q_end - this.q_start + 1;                
                if (this.q_pos_strand){
                    this.s_start += nongap_s_bases; 
                }else {
                    this.s_end -= nongap_s_bases; 
                }
                this.s_len = this.s_end - this.s_start + 1;
                
            }
        
        }catch(Exception ex){
            System.out.println(ex.getMessage());
        }
        
    }
    
    public boolean distinctQuery(BlastHSP aHsp, int overlap){
        if (!aHsp.q_contig.equals(this.q_contig))
            return true;
        if (aHsp.q_start > this.q_end  || aHsp.q_end<this.q_start)
            return true;
        
         if( aHsp.q_start<=this.q_end )
        
        if (overlap>0){
            if( (aHsp.q_start<=this.q_end && aHsp.q_start>=this.q_end-overlap) || (aHsp.q_end>=this.q_start && aHsp.q_end<=this.q_start+overlap) )
                return true;
        }
        
        
        return false;
    }
    
    public boolean addToListNoOverlap(List<BlastHSP> list, int overlap){
        // tries to add to the list if there is no major overlapping with query
        for(BlastHSP hsp: list){
            if(!distinctQuery(hsp, overlap))
                return false;            
        }
        // we get here when there is no major overlapping 
        list.add(this);
        return true;
    }

    public boolean addToList(List<BlastHSP> list, boolean log){
        //
        for(BlastHSP hsp: list){
            if (log){
            }
            boolean log2 = false;
           this.trimQuery2(hsp, log2); 
           if(log){
           }
            if (this.q_len<200)
                return false;
        }
        list.add(this);
        return true;
    }

    
    public boolean comesBefore(BlastHSP hsp, int overlap){
        if (!hsp.q_contig.equals(this.q_contig))
            return false;
        if (hsp.q_start > this.q_end)
            return true;
        if (overlap>0){
            if( hsp.q_start<=this.q_end && hsp.q_start>=this.q_end-overlap)
                return true;
        }
        
        return false;
    }
    public boolean comesAfter(BlastHSP hsp, int overlap){
        if (!hsp.q_contig.equals(this.q_contig))
            return false;
        if (hsp.q_end<this.q_start)
            return true;
        if (overlap>0){
            if( hsp.q_end>=this.q_start && hsp.q_end<=this.q_start+overlap )
                return true;
        }
        
        return false;
    }
    
    public boolean equalLength(){
        //all alignments must be equal
        if (this.mafLines==null || this.mafLines.size()==0)
            return true;
        long first_len = this.mafLines.get(0).alignment_block.length();
        for (MAFSLine line : mafLines)
            if (line.alignment_block.length()!= first_len)
                return false;
        return true;
    }
    public String toString() { 
       return String.format("%s\t%d\t%d\t%s\t%d\t%d\t%s", q_contig, q_start, q_end, s_contig, s_start,s_end, orientation==1?"+/+":"+/-");
    } 
    public void readFromLine(String line){
        //read from : contig00001	9	8080	Triticum_aestivum_Kronos_EIv1_scaffold_029728	54506	62618	98.410000
        String[] parts = line.split("\t");
        
        this.q_contig = parts[0];
        this.q_start = Integer.parseInt(parts[1]);
        this.q_end = Integer.parseInt(parts[2]);
        this.s_contig = parts[3];
        this.s_start = Integer.parseInt(parts[4]);
        this.s_end = Integer.parseInt(parts[5]);
        this.identity_percent = Double.parseDouble(parts[6])/100;
    }
    public String toBedLine(){
        //prepare subject for bed line
        int ss = 0;
        int se = 0;
        if (this.s_start>this.s_end){
            ss = this.s_end;
            se = this.s_start;
        }else{
            ss = this.s_start;
            se = this.s_end;
        }
        
        return String.format("%s\t%d\t%d", s_contig, ss-1, se);
    }
    public String toGFFLine(String chromosome){
        //prepare subject for gff line
        int ss = 0;
        int se = 0;
        if (this.s_start>this.s_end){
            ss = this.s_end;
            se = this.s_start;
            orientation = 2;
        }else{
            ss = this.s_start;
            se = this.s_end;
            orientation = 1;
        }
        
        return String.format("%s\tBlast\tHSP\t%d\t%d\t%f\t%s\t.\tQuery=%s_%d_%d_%s", s_contig, ss, se, identity_percent, orientation==1?"+":"-" , q_contig, q_start, q_end, chromosome);
    }
    
    public Set<String> getFirstContigPartSet(){
        //this is adhoc to process Cactus output
    
     Set<String> out = new HashSet();
     for(MAFSLine line : mafLines){
         String first_part = line.contig.split("\\.")[0];
         if (!first_part.startsWith("Anc0") && !first_part.startsWith("_MINIGRAPH"))
         out.add(first_part);
     }
     return out;
    }
    public Set<String> getFirstContigPartSet2(){
        //this is adhoc to process Cactus output
    
     Set<String> out = new HashSet();
     for(MAFSLine line : mafLines)
         out.add(line.contig.split("_")[0]);
     return out;
    }
    
    public void calculateOverlaps(Map<String, Integer> overlaps, Map<String, List<MAFSLine>> lists, int max_overlap){
        //this is adhoc to process Cactus output
        for(MAFSLine line : mafLines){
          String bacteria_name = line.contig.split("\\.")[0];
          int overlap = line.addToListNoOverlap(lists.get(bacteria_name), max_overlap);
          if(overlap>0){
              Integer cc = overlaps.get(bacteria_name);
              if(cc==null)
                  overlaps.put(bacteria_name, overlap);
              else
                  overlaps.put(bacteria_name, cc+overlap);
          }
              
          
          
        }
    }
    
    public void validateMAFEntry(){
        //check the size of all alignment blocks to be equal
        long len = mafLines.get(0).alignment_block.length();
        for(MAFSLine line : mafLines){
            long slen = line.alignment_block.length();
            if (slen!=len){
                System.out.println(String.format("Error: contig %s from %d to %d is not equal size with first (query) line", line.contig, line.startPOS, line.endPOS));
                System.out.println(mafLines.get(0).alignment_block);
                System.out.println(line.alignment_block);
            }
        }
    }
}

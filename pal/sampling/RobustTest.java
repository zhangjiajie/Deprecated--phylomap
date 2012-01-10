package pal.sampling;

import pal.alignment.*;
import pal.misc.Identifier;
import pal.misc.SimpleIdGroup;
import pal.tree.*;


import java.util.*;
import java.io.*;


public class RobustTest {

	/**
	 * @param args
	 */
	private ReadAlignment alnOri;
	private ReadAlignment alnExclude;
	private ReadAlignment alnCurr;
	private String output;
	private String subtree;
	private String newtree;
	
	public RobustTest (String alnOri, String alnExclude, String output){
		try{
			this.alnOri= new ReadAlignment(alnOri);
			this.alnExclude=new ReadAlignment(alnExclude);
			this.output=output;
		}catch(Exception e){
			e.printStackTrace();
		}
	}
	
	public RobustTest (String subaln, String subttree, String newtree, String output){
		try{
			this.alnExclude = new ReadAlignment(subaln);
			this.subtree=subttree;
			this.newtree=newtree;
			this.output=output;
		}catch(Exception e){
			e.printStackTrace();
		}
		
		
	}
	
	public SubAlignment excludeSequences(){
		int numExclude=alnExclude.getIdCount();
		int numOri=alnOri.getIdCount();
		int numCurr=numOri-numExclude;
		int numSites=alnOri.getSiteCount();
		
		// Reserve memory
		SimpleIdGroup idGroup = new SimpleIdGroup(numCurr);
		char[][] data = new char[numCurr][numSites];
		HashSet hs=new HashSet();
		
		for(int i=0;i<numExclude; i++){
			hs.add(alnExclude.getIdentifier(i).getName());
		}
		
		int cnt=0;
		for(int i=0; i<numOri; i++){
			if(!hs.contains(this.alnOri.getIdentifier(i).getName())){
				data[cnt]=this.alnOri.getAlignedSequenceString(i).toCharArray();
				Identifier id=alnOri.getIdentifier(i);
				idGroup.setIdentifier(cnt, id);
				cnt++;
			}
		}
		
		SubAlignment subaln=new SubAlignment(numCurr, numSites, data, idGroup);
		
		this.alnCurr=subaln;
		
		return subaln;
	}
	
	public void addSequences(int num, String alnout){
		
		int numSites=alnOri.getSiteCount();
		int numSeq=this.alnExclude.getSequenceCount()+num;
		
		SimpleIdGroup idGroup = new SimpleIdGroup(numSeq);
		char[][] data = new char[numSeq][numSites];
		
		int cnt=0;
		
		int[] perm=Permutation.next(this.alnCurr.getSequenceCount());
		
		for(int i=0; i<num; i++){
			data[cnt]=alnCurr.getAlignedSequenceString(perm[i]).toCharArray();
			//System.out.println(data[i]);
			Identifier id=alnCurr.getIdentifier(perm[i]);
			//hs.add(id.getName());
			idGroup.setIdentifier(cnt, id);
			cnt++;
		}
		
		for(int i=0;i<this.alnExclude.getIdCount();i++){
			data[cnt]=this.alnExclude.getAlignedSequenceString(i).toCharArray();
			Identifier id=alnExclude.getIdentifier(i);
			idGroup.setIdentifier(cnt, id);
			cnt++;
		}
		
		SubAlignment subaln=new SubAlignment(numSeq, numSites, data, idGroup);
		subaln.toFile(new File(alnout));
	}
	
	public void RFtest(int start, int end, int step){
		StringBuffer sb=new StringBuffer("");
		HashSet hs = new HashSet();
		for(int i=0;i<this.alnExclude.getSequenceCount();i++){
			hs.add(this.alnExclude.getIdentifier(i).getName());
		}
		
		for(int i=start; i<end; i+=step){
			ExtractSubtree es = new ExtractSubtree(this.newtree+i, hs);
			SimpleTree st=es.getSubtree();
			es.writeSubTree(new File(this.output+"newtree"+i+".tre"), st);
			double dis=TreeDistance.RFdistanceScaled(this.output+"newtree"+i+".tre", this.subtree);
			if(dis>=0){
				sb.append(i+"	"+dis);
				System.out.println(i+"	"+dis);
			}
		}
		//ExtractSubtree es = new ExtractSubtree(this.trueBigtreeFilename, hs);
		//SimpleTree st=es.getSubtree();
		//es.writeSubTree(this.trueSubtreeFile, st);
	}
	
	
	
	public void genAppendedAlignment(int start, int end, int step){
		this.excludeSequences();
		for(int i=start; i<end; i+=step){
			this.addSequences(i, this.output+""+i+".phy");
		}
	}
	
	public void RFimprove(int start, int end, int step) throws Exception {
		String folder = "/home/zhangje/storage/data/d2554/l1232/PhylomapSample/Rtest/";
		
		for(int i=start; i<end; i+=step){
			
			ReadAlignment laln=new ReadAlignment(folder + "PhylomapAln98_rt"+i+".phy");
			this.alnOri=laln;
			ReadAlignment randaln=this.excludeSequences();
			HashSet hs = new HashSet();
			for(int k=0;k<randaln.getSequenceCount();k++){
				hs.add(randaln.getIdentifier(k).getName());
			}
			
			ExtractSubtree estrue = new ExtractSubtree(folder+"TrueTree.tre", hs);
			
			ExtractSubtree esinfer = new ExtractSubtree(folder+"RAxML_bestTree.rt"+i, hs);
			
			//System.out.println (estrue.getSubtree());
			
			//System.out.println (esinfer.getSubtree());
			
			double dis=TreeUtils.getRobinsonFouldsRescaledDistance(estrue.getSubtree(), esinfer.getSubtree());
			
			if(dis>=0){
				System.out.println(i+"	"+dis);
			}
		}
	}
	
	
	public void printScript(){
		for(int i=0; i<2000; i+=10){
			System.out.println("~/bin/raxmlHPC-PTHREADS-SSE3 -T 48 -m GTRGAMMA -s PhylomapAln98_rt"+i+".phy -p 1234 -n rt"+i);
		}
	}
	
	public static void main(String[] args) throws Exception{
		// TODO Auto-generated method stub
		String aln="/home/zhangje/storage/data/d2554/l1232/RandSample/alignment_sim.phy";
		String subaln="/home/zhangje/storage/data/d2554/l1232/PhylomapSample/PhylomapAln98.phy";
		//String outfolder="/home/zhangje/storage/data/d2554/l1232/RandSample/Rtest/RandAln98_rt";
		String outfolder="/home/zhangje/storage/data/d2554/l1232/RandSample/Rtest/";
		String subtree="/home/zhangje/storage/data/d2554/l1232/RandSample/Rtest/RandTree98.tre";
		String newtree="/home/zhangje/storage/data/d2554/l1232/RandSample/Rtest/RAxML_bestTree.rt";
		RobustTest rt = new RobustTest (aln, subaln, outfolder);
		//RobustTest rt = new RobustTest (subaln, subtree, newtree, outfolder);
		//rt.genAppendedAlignment(0, 2000, 10);
		//rt.printScript();	
		//rt.RFtest(0, 850, 10);
		rt.RFimprove(10, 2000, 10);
	}

}

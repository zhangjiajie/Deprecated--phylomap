package pal.sampling;
import java.util.*;
import java.io.*;
import pal.alignment.*;
import pal.misc.*;
import pal.tree.*;

public class RandomSample {

	/**
	 * @param args
	 */
	private String trueBigtreeFilename;
	private File trueSubtreeFile;
	private File sampleAlignmentFile;
	private int numSampling=0;
	private ReadAlignment aln;
	private String raxmlmatrix;
	private String soutput;
	private String njtreeout;
	
	/*
	public RandomSample(String alignment, String bigTree, String subTree, String sampleAlignment, int numSampling){
		try{
			this.aln=new ReadAlignment(alignment);
			this.trueBigtreeFilename=bigTree;
			trueSubtreeFile = new File(subTree);
			sampleAlignmentFile = new File (sampleAlignment);
			this.numSampling=numSampling;
		}catch(Exception e){
			e.printStackTrace();
		}
	}
	
	public RandomSample(String alignment, String bigTree, String subTree, String sampleAlignment, String raxmlmatrix, String output){
		try{
			this.aln=new ReadAlignment(alignment);
			this.trueBigtreeFilename=bigTree;
			trueSubtreeFile = new File(subTree);
			sampleAlignmentFile = new File (sampleAlignment);
			//this.numSampling=numSampling;
			this.raxmlmatrix=raxmlmatrix;
			this.soutput=output;
		}catch(Exception e){
			e.printStackTrace();
		}
	}
	*/
	
	//This one is in use 
	public RandomSample(String alignment, String bigTree, String subTree, String sampleAlignment, String raxmlmatrix, String output, String njtree){
		try{
			this.aln=new ReadAlignment(alignment);
			this.trueBigtreeFilename=bigTree;
			trueSubtreeFile = new File(subTree);
			sampleAlignmentFile = new File (sampleAlignment);
			//this.numSampling=numSampling;
			this.raxmlmatrix=raxmlmatrix;
			this.soutput=output;
			this.njtreeout=njtree;
		}catch(Exception e){
			e.printStackTrace();
		}
	}
	
	/*
	public void runSampling(){
		
		int numSeqs=numSampling;
		int numSites=aln.getSiteCount();
		
		// Reserve memory
		SimpleIdGroup idGroup = new SimpleIdGroup(numSampling);
		char[][] data = new char[numSeqs][numSites];
		HashSet hs=new HashSet();
		
		int[] perm=Permutation.next(aln.getSequenceCount());
		
		for(int i=0; i<this.numSampling; i++){
			data[i]=aln.getAlignedSequenceString(perm[i]).toCharArray();
			//System.out.println(data[i]);
			Identifier id=aln.getIdentifier(perm[i]);
			hs.add(id.getName());
			idGroup.setIdentifier(i, id);
		}
		
		SubAlignment subaln=new SubAlignment(numSeqs, numSites, data, idGroup);
		subaln.toFile(this.sampleAlignmentFile);
		ExtractSubtree es = new ExtractSubtree(this.trueBigtreeFilename, hs);
		SimpleTree st=es.getSubtree();
		es.writeSubTree(this.trueSubtreeFile, st);
		
	}
	*/
	
	//This one is in use
	public void runSampling(int k){
		
		this.numSampling=k;
		int numSeqs=numSampling;
		int numSites=aln.getSiteCount();
		
		// Reserve memory
		SimpleIdGroup idGroup = new SimpleIdGroup(numSampling);
		char[][] data = new char[numSeqs][numSites];
		HashSet hs=new HashSet();
		
		int[] perm=Permutation.next(aln.getSequenceCount());
		
		for(int i=0; i<this.numSampling; i++){
			data[i]=aln.getAlignedSequenceString(perm[i]).toCharArray();
			//System.out.println(data[i]);
			Identifier id=aln.getIdentifier(perm[i]);
			hs.add(id.getName());
			idGroup.setIdentifier(i, id);
		}
		
		SubAlignment subaln=new SubAlignment(numSeqs, numSites, data, idGroup);
		subaln.toFile(this.sampleAlignmentFile);
		ExtractSubtree es = new ExtractSubtree(this.trueBigtreeFilename, hs);
		SimpleTree st=es.getSubtree();
		es.writeSubTree(this.trueSubtreeFile, st);
		
		Error err=new Error(this.raxmlmatrix, this.sampleAlignmentFile.getAbsolutePath(), this.aln);
		double errvalue=err.getAverageError();
		double pdvalue=err.getPDvalue(st);
		File foutput=new File(this.soutput);
		StringBuffer sb =new StringBuffer("");
		sb.append("samplesize	"+k+"\n");
		sb.append("error	"+errvalue+"\n");
		sb.append("pd	"+pdvalue+"\n");
		try{
			FileWriter fw=new FileWriter(foutput);
			fw.write(sb.toString());
			fw.close();
		}catch(Exception e){
			e.printStackTrace();
		}
		
		NeighborJoiningTree njt = new NeighborJoiningTree (err.getSampleDistanceMatrix());
		es.writeSubTree(new File(this.njtreeout), njt);
	}
	
	public void runSamplingOnly(int k){
		
		this.numSampling=k;
		int numSeqs=numSampling;
		int numSites=aln.getSiteCount();
		
		// Reserve memory
		SimpleIdGroup idGroup = new SimpleIdGroup(numSampling);
		char[][] data = new char[numSeqs][numSites];
		HashSet hs=new HashSet();
		
		int[] perm=Permutation.next(aln.getSequenceCount());
		
		for(int i=0; i<this.numSampling; i++){
			data[i]=aln.getAlignedSequenceString(perm[i]).toCharArray();
			//System.out.println(data[i]);
			Identifier id=aln.getIdentifier(perm[i]);
			hs.add(id.getName());
			idGroup.setIdentifier(i, id);
		}
		
		SubAlignment subaln=new SubAlignment(numSeqs, numSites, data, idGroup);
		subaln.toFile(this.sampleAlignmentFile);
		ExtractSubtree es = new ExtractSubtree(this.trueBigtreeFilename, hs);
		SimpleTree st=es.getSubtree();
		es.writeSubTree(this.trueSubtreeFile, st);
		
		/*
		Error err=new Error(this.raxmlmatrix, this.sampleAlignmentFile.getAbsolutePath(), this.aln);
		double errvalue=err.getAverageError();
		double pdvalue=err.getPDvalue(st);
		File foutput=new File(this.soutput);
		StringBuffer sb =new StringBuffer("");
		sb.append("samplesize	"+k+"\n");
		sb.append("error	"+errvalue+"\n");
		sb.append("pd	"+pdvalue+"\n");
		try{
			FileWriter fw=new FileWriter(foutput);
			fw.write(sb.toString());
			fw.close();
		}catch(Exception e){
			e.printStackTrace();
		}
		*/
		Error err=new Error(this.raxmlmatrix, this.sampleAlignmentFile.getAbsolutePath(), this.aln);
		NeighborJoiningTree njt = new NeighborJoiningTree (err.getSampleDistanceMatrix());
		es.writeSubTree(new File(this.njtreeout), njt);
		
	}
	
	
	public static void sampleManytimes(int start, int end, int step, String folder){
		/*
		 * Require 3 files in the folder:
		 * alignment_sim.phy
		 * TrueTree.tre
		 * raxmlmatrix.mx
		 */
		
		String alignment=folder+"alignment_sim.phy";
		String bigTree=folder+"TrueTree.tre";
		String subTree=folder+"RandTree";
		String sampleAlignment=folder+"RandAln";
		String raxmlmatrix=folder+"raxmlmatrix.mx";
		String output=folder+"Rand";
		String njtree=folder+"RandNJ";
		
		StringBuffer sb=new StringBuffer("");
		
		for(int i=start; i<end; i+=step){
			RandomSample rs=new RandomSample(alignment, bigTree, subTree+i+".tre", sampleAlignment+i+".phy", raxmlmatrix, output+i+".info", njtree+i+".tre");
			rs.runSampling(i);
			sb.append("~/bin/raxmlHPC-PTHREADS-SSE3 -T 48 -m GTRGAMMA -s RandAln"+i+".phy -p 1234 -n Rand"+i+"\n");
			System.out.println("Sampled: "+i);
		}
		File scriptout=new File(folder+"Rand"+start+"_"+end+".job");
		try{
			FileWriter fw=new FileWriter(scriptout);
			fw.write(sb.toString());
			fw.close();
		}catch(Exception e){
			e.printStackTrace();
		}
		
	}
	
	public static void sampleOnlyManytimes(int start, int end, int step, String folder){
		/*
		 * Require 3 files in the folder:
		 * alignment_sim.phy
		 * TrueTree.tre
		 * raxmlmatrix.mx
		 */
		
		String alignment=folder+"alignment_sim.phy";
		String bigTree=folder+"TrueTree.tre";
		String subTree=folder+"RandTree";
		String sampleAlignment=folder+"RandAln";
		String raxmlmatrix=folder+"raxmlmatrix.mx";
		String output=folder+"Rand";
		String njtree=folder+"RandNJ";
		
		StringBuffer sb=new StringBuffer("");
		
		for(int i=start; i<end; i+=step){
			RandomSample rs=new RandomSample(alignment, bigTree, subTree+i+".tre", sampleAlignment+i+".phy", raxmlmatrix, output+i+".info", njtree+i+".tre");
			rs.runSamplingOnly(i);
			sb.append("~/bin/raxmlX -m GTRGAMMA -X -s RandAln"+i+".phy -p 1234 -n Rand"+i+"\n");
			System.out.println("Sampled: "+i);
		}
		File scriptout=new File(folder+"Rand"+start+"_"+end+".job");
		try{
			FileWriter fw=new FileWriter(scriptout);
			fw.write(sb.toString());
			fw.close();
		}catch(Exception e){
			e.printStackTrace();
		}
		
	}
	
	/*
	public static void sampleManytimes(){
		String alignment="/home/zhangje/data/NP/NP_GTR_noindel_TRUE.phy";
		String bigTree="/home/zhangje/data/NP/RAxML_bestTree.np";
		String subTree="/home/zhangje/data/NP/randsample/rand";
		String sampleAlignment="/home/zhangje/data/NP/randsample/rand";
		
		StringBuffer sb=new StringBuffer("");
		
		//only to 30-1291 1501-1681
		for(int i=30; i<2000; i++){
			RandomSample rs=new RandomSample(alignment, bigTree, subTree+i+".tre", sampleAlignment+i+".phy", i);
			rs.runSampling();
			sb.append("~/bin/raxmlHPC-PTHREADS-SSE3 -T 48 -m GTRGAMMA -s ~/data/NP/randsampe/rand"+i+".phy -p 1234 -n rand"+i+"\n");
		}
		File scriptout=new File("/home/zhangje/data/NP/randsample/rand.job");
		try{
			FileWriter fw=new FileWriter(scriptout);
			fw.write(sb.toString());
			fw.close();
		}catch(Exception e){
			e.printStackTrace();
		}
		
	}
	*/
	/*
	public static void compareManytimes(){
		String trueTree="/home/zhangje/storage/data/d4114/randsample/rand";
		String inferTree="/home/zhangje/storage/data/d4114/randsample/raxml/RAxML_bestTree.rand";
		//File compareOut=new File("/home/zhangje/data/NP/randsample/30_1292.sout");
		File compareOut=new File("/home/zhangje/data/alexi/d4114/rand.sout");
		StringBuffer sb=new StringBuffer("");
		
		for(int i=30; i<3500; i++){
			double dis=TreeDistance.RFdistanceScaled(trueTree+i+".tre", inferTree+i);
			if(dis>-0.1){
				sb.append(i+"	"+dis+"\n");
				System.out.println(i+"	"+dis);
			}
		}
		
		try{
			FileWriter fw=new FileWriter(compareOut);
			fw.write(sb.toString());
			fw.close();
		}catch(Exception e){
			e.printStackTrace();
		}
		
		
		
	}
	*/
	
	public static void main(String[] args) {
		// TODO Auto-generated method stub
		RandomSample.sampleOnlyManytimes(10, 20, 10, "/lhome/lzhangje/");

	}

}

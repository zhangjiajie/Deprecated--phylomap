package pal.sampling;

import java.io.*;
import java.util.*;

import pal.alignment.*;
import pal.misc.Identifier;
import pal.misc.SimpleIdGroup;
import pal.tree.ExtractSubtree;
import pal.tree.NeighborJoiningTree;
import pal.tree.SimpleTree;

public class PhylomapSampling {

	/**
	 * @param args
	 */
	private String trueBigtreeFilename;
	private File trueSubtreeFile;
	private File sampleAlignmentFile;
	private int numSampling=0;
	private ReadAlignment aln;
	private double[][] dataMatrix;
	private String soutput;
	private String njtree;
	
	
	public PhylomapSampling(String alignment, String bigTree, String subTree, String sampleAlignment, int numSampling, double[][] dataMatrix){
		try{
			this.aln=new ReadAlignment(alignment);
			this.trueBigtreeFilename=bigTree;
			trueSubtreeFile = new File(subTree);
			sampleAlignmentFile = new File (sampleAlignment);
			this.numSampling=numSampling;
			this.dataMatrix=dataMatrix;
			//AlignmentUtils.printSequential(aln, new PrintWriter(sampleAlignment));
		}catch(Exception e){
			e.printStackTrace();
		}
	}
	
	public PhylomapSampling(String alignment, String bigTree, String subTree, String sampleAlignment, String raxmlmatrix, String output, String njtree, int numSampling){
		try{
			this.aln=new ReadAlignment(alignment);
			this.trueBigtreeFilename=bigTree;
			trueSubtreeFile = new File(subTree);
			sampleAlignmentFile = new File (sampleAlignment);
			this.numSampling=numSampling;
			this.dataMatrix=DistacneMatrix.readDistacneMatrixRaxml(new File(raxmlmatrix), alignment);
			this.soutput=output;
			this.njtree=njtree;
			//AlignmentUtils.printSequential(aln, new PrintWriter(sampleAlignment));
		}catch(Exception e){
			e.printStackTrace();
		}
	}
	
	public void runSampling(){
		
		int numSeqs=numSampling;
		int numSites=aln.getSiteCount();
		
		// Reserve memory
		SimpleIdGroup idGroup = new SimpleIdGroup(numSampling);
		char[][] data = new char[numSeqs][numSites];
		HashSet hs=new HashSet();
		
		/*
		int[] perm=Permutation.next(aln.getSequenceCount());
		
		for(int i=0; i<this.numSampling; i++){
			data[i]=aln.getAlignedSequenceString(perm[i]).toCharArray();
			//System.out.println(data[i]);
			Identifier id=aln.getIdentifier(perm[i]);
			hs.add(id.getName());
			idGroup.setIdentifier(i, id);
		}
		*/
		
		NeuralGas ng=new NeuralGas(this.dataMatrix, this.numSampling, 100);
		ng.learn();
		int codebook[]=ng.getCodebookvector();
		
		/*
		HashSet hs_codebook=new HashSet();
		for(int i=0; i<codebook.length; i++){
			hs_codebook.add(new Integer(codebook[i]));
		}
		*/
		
		for(int i=0; i<this.numSampling; i++){
			System.out.print(codebook[i]+",");
			data[i]=aln.getAlignedSequenceString(codebook[i]).toCharArray();
			Identifier id=aln.getIdentifier(codebook[i]);
			hs.add(id.getName());
			idGroup.setIdentifier(i, id);
		}
		
		SubAlignment subaln=new SubAlignment(numSeqs, numSites, data, idGroup);
		subaln.toFile(this.sampleAlignmentFile);
		ExtractSubtree es = new ExtractSubtree(this.trueBigtreeFilename, hs);
		SimpleTree st=es.getSubtree();
		es.writeSubTree(this.trueSubtreeFile, st);
		
		Error err=new Error(this.dataMatrix, this.sampleAlignmentFile.getAbsolutePath(), this.aln);
		double errvalue=err.getAverageError();
		double pdvalue=err.getPDvalue(st);
		File foutput=new File(this.soutput);
		StringBuffer sb =new StringBuffer("");
		sb.append("samplesize	"+this.numSampling+"\n");
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
		es.writeSubTree(new File(this.njtree), njt);
	}
	
	public static void sampleManytimes4114() throws Exception {
		String alignment="/home/zhangje/data/alexi/d4114/4114_GTR_sim_TRUE.phy";
		String bigTree="/home/zhangje/data/alexi/d4114/RAxML_bestTree.4114";
		String subTree="/home/zhangje/data/alexi/d4114/phylomapsample/phylomap";
		String sampleAlignment="/home/zhangje/data/alexi/d4114/phylomapsample/phylomap";
		
		String sdism="/home/zhangje/data/alexi/d4114/RAxML_distances.dism";
		double[][] data= DistacneMatrix.readDistacneMatrixRaxml(new File(sdism),alignment);
		StringBuffer sb=new StringBuffer("");
		
		//only to 30-1291 1501-1681
		for(int i=30; i<200; i++){
			PhylomapSampling rs=new PhylomapSampling(alignment, bigTree, subTree+i+".tre", sampleAlignment+i+".phy", i, data);
			rs.runSampling();
			sb.append("~/bin/raxmlHPC-PTHREADS-SSE3 -T 48 -m GTRGAMMA -s /home/zhangje/data/alexi/d4114/randsample/rand"+i+".phy -p 1234 -n rand"+i+"\n");
		}
		
		/*
		File scriptout=new File("/home/zhangje/data/alexi/d4114/randsample/rand.job");
		try{
			FileWriter fw=new FileWriter(scriptout);
			fw.write(sb.toString());
			fw.close();
		}catch(Exception e){
			e.printStackTrace();
		}
		*/

	}
	
	public static void genRunscript(){
		
		StringBuffer sb = new StringBuffer("");
		for(int i=30; i<1003; i++){
			File ftest=new File("/home/zhangje/storage/data/d4114/phylomapsample/phylomap"+i+".phy");
			if (ftest.exists()){
				sb.append("~/bin/raxmlHPC-SSE3 -m GTRGAMMA -s /hits/sco/zhangje/data/d4114/phylomapsample/phylomap"+i+".phy -p 1234 -n phylomap"+i+" \n");
			}
		}
		
		
		File scriptout=new File("/home/zhangje/data/alexi/d4114/phylomap_4114.job");
		try{
			FileWriter fw=new FileWriter(scriptout);
			fw.write(sb.toString());
			fw.close();
		}catch(Exception e){
			e.printStackTrace();
		}
		
	}
	
	
	public static void compareManytimes(){
		/*
		String trueTree="/home/zhangje/storage/data/d4114/phylomapsample/phylomap";
		String inferTree="/home/zhangje/storage/data/d4114/phylomapsample/job/RAxML_bestTree.phylomap";
		//File compareOut=new File("/home/zhangje/data/NP/randsample/30_1292.sout");
		File compareOut=new File("/home/zhangje/data/alexi/d4114/phylomap.sout");
		*/
		
		
		String trueTree="/home/zhangje/storage/data/NP/phylomapsample/phylomap";
		String inferTree="/home/zhangje/storage/data/NP/phylomapsample/job/RAxML_bestTree.phylomap";
		//File compareOut=new File("/home/zhangje/data/NP/randsample/30_1292.sout");
		File compareOut=new File("/home/zhangje/data/NP/phylomap.sout");
		
		
		
		StringBuffer sb=new StringBuffer("");
		
		for(int i=551; i<556; i++){
			double dis=TreeDistance.RFdistanceScaled(trueTree+i+".tre", inferTree+i);
			if(dis>-0.1){
				sb.append(i+"	"+dis+"\n");
				System.out.println(i+"	"+dis);
			}
		}
		
		try{
			//FileWriter fw=new FileWriter(compareOut);
			//fw.write(sb.toString());
			//fw.close();
		}catch(Exception e){
			e.printStackTrace();
		}
		
		
		
	}
	
	public static void main(String[] args) throws Exception{
		// TODO Auto-generated method stub
		String alignment="/home/zhangje/data/alexi/d4114/4114_GTR_sim_TRUE.phy";
		String bigTree="/home/zhangje/data/alexi/d4114/RAxML_bestTree.4114";
		String subTree="/home/zhangje/data/alexi/d4114/phylomapsample/sample.tre";
		String sampleAlignment="/home/zhangje/data/alexi/d4114/phylomapsample/sample.phy";
		File fdism=new File("/home/zhangje/data/alexi/d4114/phylomapsample/outfile");
		//double[][] data= DistacneMatrix.readDistanceMatrix(fdism);
		//String sdism="/home/zhangje/data/alexi/d4114/RAxML_distances.dism";
		//double[][] data= DistacneMatrix.readDistacneMatrixRaxml(new File(sdism),alignment);
		
		//PhylomapSampling ps=new PhylomapSampling(alignment, bigTree, subTree, sampleAlignment, 150, data);
		//ps.runSampling();
		//PhylomapSampling.sampleManytimes4114();
		//PhylomapSampling.genRunscript();
		//PhylomapSampling.compareManytimes();
		
		
		
	}

}

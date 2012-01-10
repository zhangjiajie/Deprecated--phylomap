package pal.sampling;
import java.io.*;
import java.util.*;

import pal.alignment.ReadAlignment;
import pal.alignment.SubAlignment;
import pal.misc.Identifier;
import pal.misc.SimpleIdGroup;
import pal.tree.ExtractSubtree;
import pal.tree.NeighborJoiningTree;
import pal.tree.SimpleTree;



public class SeqVecSampling {

	/**
	 * @param args
	 */
	double[] A= {0.3535533905932738, 0.5, 0.0};
	double[] G= {0.3535533905932738, -0.5, 0.0};
	double[] T=	{-0.3535533905932738, 0.0, 0.5};
	double[] U=	{-0.3535533905932738, 0.0, 0.5};
	double[] C=	{-0.3535533905932738, 0.0, -0.5};
	double[] Space= {0.0, 0.0, 0.0};
	HashMap<String, double[]> hm = new HashMap<String, double[]>();
	
	private String trueBigtreeFilename;
	private File trueSubtreeFile;
	private File sampleAlignmentFile;
	private int numSampling=0;
	private ReadAlignment aln;
	private double[][] dataMatrix;
	private String raxmlmatrix;
	private String soutput;
	private String njtree;
	
	/*
	public SeqVecSampling(String alignment, String bigTree, String subTree, String sampleAlignment, int numSampling){
		hm.put("A", A);
		hm.put("a", A);
		hm.put("G", G);
		hm.put("g", G);
		hm.put("T", T);
		hm.put("t", T);
		hm.put("U", U);
		hm.put("u", U);
		hm.put("C", C);
		hm.put("c", C);
		hm.put("-", Space);
		try{
			this.aln=new ReadAlignment(alignment);
			this.trueBigtreeFilename=bigTree;
			trueSubtreeFile = new File(subTree);
			sampleAlignmentFile = new File (sampleAlignment);
			this.numSampling=numSampling;
			this.genDismatrix();
		}catch(Exception e){
			e.printStackTrace();
		}
	}
	*/
	
	public SeqVecSampling(String alignment, String bigTree, String subTree, String sampleAlignment, String raxmlmatrix, String output, String njtree, int numSampling){
		hm.put("A", A);
		hm.put("a", A);
		hm.put("G", G);
		hm.put("g", G);
		hm.put("T", T);
		hm.put("t", T);
		hm.put("U", U);
		hm.put("u", U);
		hm.put("C", C);
		hm.put("c", C);
		hm.put("-", Space);
		try{
			this.aln=new ReadAlignment(alignment);
			this.trueBigtreeFilename=bigTree;
			trueSubtreeFile = new File(subTree);
			sampleAlignmentFile = new File (sampleAlignment);
			this.numSampling=numSampling;
			this.raxmlmatrix=raxmlmatrix;
			this.soutput=output;
			this.genDismatrix();
			this.njtree=njtree;
		}catch(Exception e){
			e.printStackTrace();
		}
	}
	
	public double[] translate(String seq){
		int len=seq.trim().length()*3;
		double[] vec=new double[len];
		int size = seq.trim().length();
		String nn;
		int cnt=0;
		double[] v;
		for(int i=0; i<size; i++){
			nn=seq.substring(i,i+1);
			v= hm.get(nn);
			if(v==null){
				v=Space;
			}
			for(int k=0;k<v.length;k++){
				vec[cnt]=v[k];
				cnt++;
			}
		}
		return vec;
	}
	
	public double[][] genDismatrix (){
		
			int numSeq=aln.getSequenceCount();
			int numSite=aln.getSiteCount();
			dataMatrix = new double[numSeq][numSite*3];
			for(int i=0; i<numSeq; i++){
				String seq=aln.getAlignedSequenceString(i);
				double[] vec=this.translate(seq);
				dataMatrix[i]=vec;
				//System.out.println("Translating sequence: "+i);
			}
			return dataMatrix;

	}

	public void runSampling(){
		
		int numSeqs=numSampling;
		int numSites=aln.getSiteCount();
		
		// Reserve memory
		SimpleIdGroup idGroup = new SimpleIdGroup(numSampling);
		char[][] data = new char[numSeqs][numSites];
		HashSet hs=new HashSet();
		
		NeuralGas ng=new NeuralGas(this.dataMatrix, this.numSampling, 100);
		ng.learn();
		int codebook[]=ng.getCodebookvector();
		
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
		
		Error err=new Error(this.raxmlmatrix, this.sampleAlignmentFile.getAbsolutePath(), this.aln);
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
	
	public static void main(String[] args) {
		// TODO Auto-generated method stub
		String alignment="/home/zhangje/storage/data/NP/NP_GTR_noindel_TRUE.phy";
		String bigTree="/home/zhangje/storage/data/NP/bigTree.tre";
		String subTree="/home/zhangje/storage/data/NP/svsample/sample.tre";
		String sampleAlignment="/home/zhangje/storage/data/NP/svsample/sample.phy";
		//SeqVecSampling svs = new SeqVecSampling(alignment, bigTree, subTree, sampleAlignment, 30);
		//svs.runSampling();
	}

}

package pal.sampling;

import java.io.*;
import java.util.*;

import pal.alignment.*;
import pal.misc.Identifier;
import pal.misc.SimpleIdGroup;
import pal.tree.ExtractSubtree;
import pal.tree.NeighborJoiningTree;
import pal.tree.SimpleTree;

public class PamSampling {

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
	
	public PamSampling(String alignment, String bigTree, String subTree, String sampleAlignment, String raxmlmatrix, String output, String njtree, int numSampling){
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
	
	public double runSampling(){
		
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
		
		//NeuralGas ng=new NeuralGas(this.dataMatrix, this.numSampling, 100);
		Date d1=new Date();
		PAM pam=new PAM(this.dataMatrix, this.numSampling);
		//ng.learn();
		int codebook[]=pam.getCodebook();
		Date d2=new Date();
		
		double lasttime=((double)(d2.getTime()-d1.getTime()))/60000.0;
		
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
		sb.append("running time	"+lasttime+" minutes \n");
		try{
			FileWriter fw=new FileWriter(foutput);
			fw.write(sb.toString());
			fw.close();
		}catch(Exception e){
			e.printStackTrace();
		}
		
		NeighborJoiningTree njt = new NeighborJoiningTree (err.getSampleDistanceMatrix());
		es.writeSubTree(new File(this.njtree), njt);
		return lasttime;
	}
	
	
	public static void main(String[] args) {
		// TODO Auto-generated method stub
		String alignment="/lhome/lzhangje/storage/data/d2554/l1232/alignment_sim.phy";
		String bigTree="/lhome/lzhangje/storage/data/d2554/l1232/TrueTree.tre";
		String subTree="/lhome/lzhangje/pam10.subtree";
		String sampleAlignment="/lhome/lzhangje/pam10.subalign";
		String raxmlmatrix="/lhome/lzhangje/storage/data/d2554/l1232/raxmlmatrix.mx";
		String output="/lhome/lzhangje/pam10.out";
		String njtree="/lhome/lzhangje/pam10.njtree";
		
		
		PamSampling ps = new PamSampling(alignment, bigTree, subTree, sampleAlignment, raxmlmatrix, output, njtree, 10);
		double lasttime=ps.runSampling();
		System.out.println();
		System.out.println("Clustering running for: " + lasttime +" minutes");
		
	}

}

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

public class CVSampling {

	/**
	 * @param args
	 */
	private String trueBigtreeFilename;
	private File trueSubtreeFile;
	private File sampleAlignmentFile;
	private int numSampling=0;
	private ReadAlignment aln;
	private double[][] dataMatrix;
	private HashSet hs =new HashSet();
	private String[] kmer;
	private int k=5;
	private String raxmlmatrix;
	private String soutput;
	private String njtree;
	
	/*
	public CVSampling(String alignment, String bigTree, String subTree, String sampleAlignment, int numSampling, int k){

		try{
			this.aln=new ReadAlignment(alignment);
			this.trueBigtreeFilename=bigTree;
			trueSubtreeFile = new File(subTree);
			sampleAlignmentFile = new File (sampleAlignment);
			this.numSampling=numSampling;
			this.genVectorMatrix();
			this.k=k;
		}catch(Exception e){
			e.printStackTrace();
		}
	}
	*/
	
	public CVSampling(String alignment, String bigTree, String subTree, String sampleAlignment, String raxmlmatrix, String output, String njtree, int numSampling, int k){

		try{
			this.aln=new ReadAlignment(alignment);
			this.trueBigtreeFilename=bigTree;
			trueSubtreeFile = new File(subTree);
			sampleAlignmentFile = new File (sampleAlignment);
			this.numSampling=numSampling;
			this.genVectorMatrix();
			this.k=k;
			this.raxmlmatrix=raxmlmatrix;
			this.soutput=output;
			this.njtree=njtree;
		}catch(Exception e){
			e.printStackTrace();
		}
	}
	
	public static void convertToCVtreeFormatc (File fin){
		try{
			FileWriter fw=new FileWriter(new File(fin.getAbsolutePath()+".cv"));
			FileReader fr=new FileReader(fin);
			BufferedReader bf=new BufferedReader(fr);
			String s=bf.readLine();
			while(s!=null){
				if(s.startsWith(">")){
					s=s.trim();
					String added = s.substring(1);
					fw.write(s+"|"+added+"\n");
				}else{
					fw.write(s.trim()+"\n");
				}
				s=bf.readLine();
			}
			
			
			bf.close();
			fr.close();
			fw.close();
			
			
		}catch(Exception e){
			e.printStackTrace();
		}
		
		
	}
	
	//scan one sequence for possible kwords
	private HashSet<String> scanOneSeq(HashSet<String> hs, String seq, int len){
		seq=new String (seq.replaceAll("-", ""));
		for(int i=0;i<seq.length()-len+1;i++){
			hs.add(seq.substring(i, i+len));
		}
		return hs;
	}
	
	private HashSet<String> scanSeqSet(HashSet<String> hs, List<String> seqlist, int len){
		for(int i=0; i<seqlist.size(); i++){
			hs=this.scanOneSeq(hs, seqlist.get(i), len);
		}
		return hs;
	}
	
	public double[] seq2vec (String seq, String[] kmer, int k){
		seq=new String (seq.replaceAll("-", ""));
		int size = kmer.length;
		int len = seq.length();
		int denom = len-k+1;
		double[] vec = new double[size];
		String[] skmer=new String[denom];
		for(int i=0; i<denom; i++){
			skmer[i]=seq.substring(i, i+k);
		}
		String currkmer;
		for(int i=0; i<size; i++){
			double cnt=0;
			currkmer=kmer[i];
			for(int j=0;j<denom;j++){
				if(currkmer.equals(skmer[j])){
					cnt+=1.0;
				}
			}
			vec[i]=cnt/((double)(denom));
		}
		return vec;
	}
	
	
	private double[][] genVectorMatrix(){
		int numSeq=aln.getSequenceCount();
		int numSite=aln.getSiteCount();
		for(int i=0; i<numSeq; i++){
			String seq=aln.getAlignedSequenceString(i);
			this.hs=this.scanOneSeq(this.hs, seq, this.k);
		}
		Object[] objs= this.hs.toArray();
		this.kmer=new String[objs.length];
		for(int i=0; i<hs.size(); i++){
			this.kmer[i]=(String)objs[i];
		}
	
		//this.kmer=(String[])this.hs.toArray();
		int size=this.kmer.length;
		this.dataMatrix=new double[numSeq][size];
		for(int i=0; i<numSeq; i++){
			String seq=aln.getAlignedSequenceString(i);
			this.dataMatrix[i]=this.seq2vec(seq, this.kmer, this.k);
		}
		
		return this.dataMatrix;
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
		//File fin = new File("/home/zhangje/data/NP/cvsample/NP_sim1.fas");
		
		
		String alignment="/home/zhangje/storage/data/NP/NP_GTR_noindel_TRUE.phy";
		String bigTree="/home/zhangje/storage/data/NP/bigTree.tre";
		String subTree="/home/zhangje/storage/data/NP/cvsample/sample.tre";
		String sampleAlignment="/home/zhangje/storage/data/NP/cvsample/sample.phy";
		//CVSampling svs = new CVSampling(alignment, bigTree, subTree, sampleAlignment, 30, 10);
		//svs.runSampling();
		String s="t2556     ATTTT-TATTCAACC---TTGACTTGGTGCGCT-TGGAGCGAAATCC--GCGTTAAAATTAGCATCCACTATTG-ATT--------------------------------------------GTA-ATCCGACA-AGATTATTACA--GGAAGACGAA-TTTACATTGATTTGTTGATAGATACGTG-ATATATTACTGGCG--------T-GGATTTGGGCTTCAAGACGTGC--TATCATCTCG---AGATAG----------------------TTTACT-TGCTCAGTCTCTAGT--GAGTTGTA----AAGGG---------A-AG-CTGCTTG-TCAAGTTGTCT------CTGCC-----GTGGATC--TAATTGGAA-TT-AA-A--CATTTTTGTC-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------TGT--C--A-GT-ACCAGGTTATAAGATTCCAA-GATCGTG---ACT---AGGACG--------GGCCAGTAT--A--GCGTTAACGG----------------------------------------------------------------------------------------------------G-ACGGTT---G---TTC-CCTTGT-CCAT-TTG-ACTT--CTT--CCCAGCACGTGTCAT----GGGGGTCAAGCACATCCAATATATCCTTCCGCAGT------------AATG-ATG---ACTTGACTACCTCTTGAATAGC-AGCGGACC-TCGGA----ACTAGGGCGTC-CAAGTCCCG---TAATCT--TAA----GACATCACAA--CGCGCAGATGGC-----TCGGG------------------------TGC-----TATCAATATTCGAGACTTTAGCGAGGACCCGAGAT-TAT-AGAGATTCTTGTGATGTAGTATTA----CTTC--TTTC--AC-TGACATCTCTATAAGCGC-ATG--GATTGTTCTCAGTTC-CC---TTTAGTTAAATGAATGT-AGATC-TT-TGGAA-------GTAAACAGCTCGGT-------------GATCTGGAATAGTGTGCGTGCC---ATT----------TCTG--GC------------------TAAGCTTAAAAGG--A--TGT------TTGGACTTGTACT-C---CTAGAGA--CTA--AT-TCTAT-------------T--TGGGGGGTGT-GT--CCAGCGTG-TGCGAACGTTTC-CTTACGA-TC--GTGGTATGAGATGG------------CCGGAGACAGTGCAAAT--GCGGGTTCGGGAT------------------GCAATGTCTGT-TCAACTTGAGGTCATAGT-C--ATTAC--CTTTATTAAA-GCTTGGAGCTTAGAGGATTAAGGTGTACTCGT-----------CATGAATACTACTGGACTTCATGGAGGCGTATGTATT-TA-T-CATACGCGCAAATTAAAAA-TA----------------------------TA--TTAATTCCGT--TGATGTGCA-----------------------AT---GT     ";
		String s2=new String(s.replaceAll("-", ""));
		System.out.println(s);
		System.out.println(s2);
	}

}

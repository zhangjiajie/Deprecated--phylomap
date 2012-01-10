package pal.sampling;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.util.HashSet;
import java.util.List;

import pal.alignment.ReadAlignment;
import pal.alignment.SubAlignment;
import pal.misc.Identifier;
import pal.misc.SimpleIdGroup;
import pal.tree.ExtractSubtree;
import pal.tree.SimpleTree;

public class CVmSampling {

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
	private int mismatich=1;
	
	public CVmSampling(String alignment, String bigTree, String subTree, String sampleAlignment, int numSampling, int k, int mismatch){

		try{
			this.aln=new ReadAlignment(alignment);
			this.trueBigtreeFilename=bigTree;
			trueSubtreeFile = new File(subTree);
			sampleAlignmentFile = new File (sampleAlignment);
			this.numSampling=numSampling;
			this.genDistanceMatrix();
			this.k=k;
			this.mismatich=mismatch;
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
				if(CVmSampling.misMatchOne(currkmer, skmer[j])){
					//if(currkmer.equals(skmer[j])){	
					cnt+=1.0;
				}
			}
			vec[i]=cnt/((double)(denom));
		}
		return vec;
	}
	
	private double[][] genDistanceMatrix(){
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
		
	}
	
	public static boolean misMatchOne(String s, String t){
		boolean flag=false;
		for(int i=0;i<s.length();i++){
			String regx=s.substring(0,i)+"."+s.substring(i+1);
			if(t.matches(regx)){
				flag = true;
				break;
			}
		}
		return flag; 
	}
	
	public static void main(String[] args) {
		String s="abcdef";
		//System.out.println(s.matches("ab.ef"));
		//System.out.println(CVmSampling.misMatchOne(s, "abc-ef"));
		String alignment="/home/zhangje/storage/data/NP/NP_GTR_noindel_TRUE.phy";
		String bigTree="/home/zhangje/storage/data/NP/bigTree.tre";
		String subTree="/home/zhangje/storage/data/NP/cvsample/sample.tre";
		String sampleAlignment="/home/zhangje/storage/data/NP/cvsample/sample.phy";
		CVmSampling svs = new CVmSampling(alignment, bigTree, subTree, sampleAlignment, 30, 9, 1);
		svs.runSampling();
		
	}
}

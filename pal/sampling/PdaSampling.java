package pal.sampling;
import java.io.*;
import java.util.*;
import pal.alignment.*;
import pal.misc.Identifier;
import pal.misc.SimpleIdGroup;
import pal.tree.*;

public class PdaSampling {

	/**
	 * @param args
	 */
	public static void extracAlignment(String saln, String staxa, String sout){
		try{
			int cnt=0;
			ReadAlignment aln = new ReadAlignment(saln);
			HashSet<String> hs=new HashSet<String>();
			//FileWriter fw=new FileWriter(new File(sout));
			FileReader fr=new FileReader(staxa);
			BufferedReader bf=new BufferedReader(fr);
			String s=bf.readLine();
			while(s!=null){
				cnt++;
				s=s.trim();
				hs.add(s);
				s=bf.readLine();
			}
			bf.close();
			fr.close();
			int numSeqs=cnt;
			int numSites=aln.getSiteCount();
			
			// Reserve memory
			SimpleIdGroup idGroup = new SimpleIdGroup(numSeqs);
			char[][] data = new char[numSeqs][numSites];
			int idx=0;
			
			//screen the whole data set 
			for(int i=0; i<aln.getSequenceCount(); i++){
				if(hs.contains(aln.getIdentifier(i).getName())){
					data[idx]=aln.getAlignedSequenceString(i).toCharArray();
					Identifier id=aln.getIdentifier(i);
					idGroup.setIdentifier(idx, id);
					idx++;
				}
			}
			//fw.write(sb.toString());
			//fw.close();
			
			SubAlignment subaln=new SubAlignment(numSeqs, numSites, data, idGroup);
			subaln.toFile(new File(sout));
			
			
		}catch(Exception e){
			e.printStackTrace();
		}
	}
	
	public static void extracAlignment(String saln, String staxa, String stree, String sout, String raxmlmatrix, String output, String njtree, int numSampling){
		try{
			int cnt=0;
			ReadAlignment aln = new ReadAlignment(saln);
			HashSet<String> hs=new HashSet<String>();
			//FileWriter fw=new FileWriter(new File(sout));
			FileReader fr=new FileReader(staxa);
			BufferedReader bf=new BufferedReader(fr);
			String s=bf.readLine();
			while(s!=null){
				cnt++;
				s=s.trim();
				hs.add(s);
				s=bf.readLine();
			}
			bf.close();
			fr.close();
			int numSeqs=cnt;
			int numSites=aln.getSiteCount();
			
			// Reserve memory
			SimpleIdGroup idGroup = new SimpleIdGroup(numSeqs);
			char[][] data = new char[numSeqs][numSites];
			int idx=0;
			
			//screen the whole data set 
			for(int i=0; i<aln.getSequenceCount(); i++){
				if(hs.contains(aln.getIdentifier(i).getName())){
					data[idx]=aln.getAlignedSequenceString(i).toCharArray();
					Identifier id=aln.getIdentifier(i);
					idGroup.setIdentifier(idx, id);
					idx++;
				}
			}
			//fw.write(sb.toString());
			//fw.close();
			
			SubAlignment subaln=new SubAlignment(numSeqs, numSites, data, idGroup);
			subaln.toFile(new File(sout));
			
			Error err=new Error(raxmlmatrix, sout, aln);
			double errvalue=err.getAverageError();
			SimpleTree st =new ReadTree(stree);
			double pdvalue=err.getPDvalue(st);
			File foutput=new File(output);
			StringBuffer sb =new StringBuffer("");
			sb.append("samplesize	"+numSampling+"\n");
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
			ExtractSubtree.writeSubTree(new File(njtree), njt);
			
			
		}catch(Exception e){
			e.printStackTrace();
		}
	}
	
	public static void genScriptPda(){
		//~/bin/pda bigTree.tre -k 29
		
		StringBuffer sb = new StringBuffer("");
		for(int i=30; i<500; i++){
			/*
			File ftest=new File("/home/zhangje/storage/data/d4114/phylomapsample/phylomap"+i+".phy");
			if (ftest.exists()){
				sb.append("~/bin/raxmlHPC-SSE3 -m GTRGAMMA -s /hits/sco/zhangje/data/d4114/phylomapsample/phylomap"+i+".phy -p 1234 -n phylomap"+i+" \n");
			}
			*/
			sb.append("~/bin/pda bigTree.tre -k "+i+"\n");
			
		}
		
		
		File scriptout=new File("/home/zhangje/storage/data/NP/pdasample/pda.script");
		try{
			FileWriter fw=new FileWriter(scriptout);
			fw.write(sb.toString());
			fw.close();
		}catch(Exception e){
			e.printStackTrace();
		}
	}
	

	
	public static void batchExtractAlignment(String folder, int start, int end, int step){
		String staxa=folder+"TrueTree.tre."; //.100.pdtaxa";
		String saln=folder+"alignment_sim.phy";
		String sout=folder+"pdaAln";
		String stree=folder+"TrueTree.tre.";
		String output=folder+"pda";
		String raxmlmatrix=folder+"raxmlmatrix.mx";
		String njtree=folder+"pdaNJ";
		for(int i=start;i<end; i+=step){
			PdaSampling.extracAlignment(saln, staxa+i+".pdtaxa", stree+i+".pdtree",sout+i+".phy", raxmlmatrix, output+i+".info",njtree+i+".tre", i);
			System.out.println("Sampled: "+i);
		}
		
	}
	
	public static void compareManytimes(){
		
		String trueTree="/home/zhangje/storage/data/d4114/pdasample/bigTree.tre.";
		String inferTree="/home/zhangje/storage/data/d4114/pdasample/job/RAxML_bestTree.pda.";
		File compareOut=new File("/home/zhangje/data/alexi/d4114/pda.sout");
	
		StringBuffer sb=new StringBuffer("");
		
		for(int i=30; i<400; i++){
			double dis=TreeDistance.RFdistanceScaled(trueTree+i+".pdtree", inferTree+i);
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
	
	public static void main(String[] args) {
		// TODO Auto-generated method stub
		String staxa="/home/zhangje/data/alexi/d4114/RAxML_bestTree.4114.100.pdtaxa";
		String saln="/home/zhangje/data/alexi/d4114/4114_GTR_sim_TRUE.phy";
		String sout="/home/zhangje/data/alexi/d4114/pdatest.txt";
		//PdaSampling.extracAlignment(saln, staxa, sout);
		//PdaSampling.genScriptPda();
		//PdaSampling.batchExtractAlignment(300, 400);
		PdaSampling.batchExtractAlignment("/home/zhangje/storage/data/d2554/l2000/PDSample/", 200, 2000, 10);
	}

}

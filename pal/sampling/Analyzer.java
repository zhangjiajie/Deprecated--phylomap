package pal.sampling;
import java.util.*;
import java.io.*;
import pal.tree.*;

public class Analyzer {

	/**
	 * @param args
	 */
	private int start=1;
	private int end=2000;
	public void summaryOutputFiles (String folder){
	//search through the folder	
		try{
			
			StringBuffer sbe=new StringBuffer("");
			StringBuffer sbpd=new StringBuffer("");
			
			for(int i= this.start; i<this.end; i++){
				File fin=new File(folder+i+".info");
				if(fin.exists()){
					FileReader fr=new FileReader(fin);
					BufferedReader bf=new BufferedReader(fr);
					String s1=bf.readLine();
					String s2=bf.readLine();
					String s3=bf.readLine();
					bf.close();
					fr.close();
					String size=s1.split("\\s")[1];
					String error=s2.split("\\s")[1];
					String pd=s3.split("\\s")[1];
					sbe.append(size+"	"+error+"\n");
					sbpd.append(size+"	"+pd+"\n");
					System.out.println(s1);
				}
			}
			
			FileWriter fw=new FileWriter(folder+".err");
			fw.write(sbe.toString());
			fw.close();
			
			FileWriter fw2=new FileWriter(folder+".pd");
			fw2.write(sbpd.toString());
			fw2.close();
		}catch(Exception e){
			e.printStackTrace();
		}
	}
	
	public void summaryRF(String folder, String stype){
		String trueTree=folder+stype+"Tree";
		String inferTree=folder+"RAxML_bestTree."+stype;
		//File compareOut=new File("/home/zhangje/data/NP/randsample/30_1292.sout");
		File compareOut=new File(folder+stype+".rf");
		StringBuffer sb=new StringBuffer("");
		
		for(int i=start; i<end; i++){
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
	
	public void summaryRFX(String folder, String stype){
		String trueTree=folder+"Rand"+"Tree";
		String inferTree=folder+"RAxML_bestTree."+stype;
		//File compareOut=new File("/home/zhangje/data/NP/randsample/30_1292.sout");
		File compareOut=new File(folder+stype+".rf");
		StringBuffer sb=new StringBuffer("");
		
		for(int i=start; i<end; i++){
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
	
	public void summaryRF_PDA(String folder){
		String trueTree=folder+"TrueTree.tre.";
		String inferTree=folder+"RAxML_bestTree.pda";
		//File compareOut=new File("/home/zhangje/data/NP/randsample/30_1292.sout");
		File compareOut=new File(folder+"pda.rf");
		StringBuffer sb=new StringBuffer("");
		
		for(int i=start; i<end; i++){
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
	
	public void summaryRF_NJ(String folder, String stype){
		String trueTree=folder+stype+"Tree";
		String inferTree=folder+stype+"NJ";
		//File compareOut=new File("/home/zhangje/data/NP/randsample/30_1292.sout");
		File compareOut=new File(folder+stype+".NJ.rf");
		StringBuffer sb=new StringBuffer("");
		
		for(int i=start; i<end; i++){
			double dis=TreeDistance.RFdistanceScaled(trueTree+i+".tre", inferTree+i+".tre");
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
	
	
	public void summaryRF_NJ_PDA(String folder){
		String trueTree=folder+"TrueTree.tre.";
		String inferTree=folder+"pda"+"NJ";
		//File compareOut=new File("/home/zhangje/data/NP/randsample/30_1292.sout");
		File compareOut=new File(folder+"pda.NJ.rf");
		StringBuffer sb=new StringBuffer("");
		
		for(int i=start; i<end; i++){
			double dis=TreeDistance.RFdistanceScaled(trueTree+i+".pdtree", inferTree+i+".tre");
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
	
	
	public void summaryRF_Par(String folder, String stype){
		
		String trueTree=folder+stype+"Tree";
		String inferTree=folder+"RAxML_parsimonyTree."+stype+"P";
		File compareOut=new File(folder+stype+".Par.rf");
		StringBuffer sb=new StringBuffer("");
		if("pda".equals(stype)){
			trueTree=folder+"TrueTree.tre.";
			for(int i=start; i<end; i++){
				double dis=TreeDistance.RFdistanceScaled(trueTree+i+".pdtree", inferTree+i+".0");
				if(dis>-0.1){
					sb.append(i+"	"+dis+"\n");
					System.out.println(i+"	"+dis);
				}
			}
		}else{
			for(int i=start; i<end; i++){
				double dis=TreeDistance.RFdistanceScaled(trueTree+i+".tre", inferTree+i+".0");
				if(dis>-0.1){
					sb.append(i+"	"+dis+"\n");
					System.out.println(i+"	"+dis);
				}
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
	
	public void summarySkewness(String folder, String stype) throws Exception{
		String trueTree=folder+stype+"Tree";
		File compareOut=new File(folder+stype+".skewness");
		StringBuffer sb=new StringBuffer("");
		if("pda".equals(stype)){
			trueTree=folder+"TrueTree.tre.";
			for(int i=start; i<end; i++){
				File f = new File(trueTree+i+".pdtree");
				if(f.exists()){
					ReadTree rtree= new ReadTree (trueTree+i+".pdtree");
					double sk=Stat.treeSkewness(rtree);
					sb.append(i+"	"+sk+"\n");
				}
			}
		}else{
			for(int i=start; i<end; i++){
				File f = new File(trueTree+i+".tre");
				if(f.exists()){
					ReadTree rtree= new ReadTree (trueTree+i+".tre");
					double sk=Stat.treeSkewness(rtree);
					sb.append(i+"	"+sk+"\n");
				}
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
	
	public void summaryVar(String folder, String stype) throws Exception{
		String trueTree=folder+stype+"Tree";
		File compareOut=new File(folder+stype+".var");
		StringBuffer sb=new StringBuffer("");
		if("pda".equals(stype)){
			trueTree=folder+"TrueTree.tre.";
			for(int i=start; i<end; i++){
				File f = new File(trueTree+i+".pdtree");
				if(f.exists()){
					ReadTree rtree= new ReadTree (trueTree+i+".pdtree");
					double sk=Stat.treeVariance(rtree);
					sb.append(i+"	"+sk+"\n");
				}
			}
		}else{
			for(int i=start; i<end; i++){
				File f = new File(trueTree+i+".tre");
				if(f.exists()){
					ReadTree rtree= new ReadTree (trueTree+i+".tre");
					double sk=Stat.treeVariance(rtree);
					sb.append(i+"	"+sk+"\n");
				}
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
	
	public void summaryMinBranch(String folder, String stype) throws Exception{
		String trueTree=folder+stype+"Tree";
		File compareOut=new File("/lhome/lzhangje/"+stype+".min");
		StringBuffer sb=new StringBuffer("");
		if("pda".equals(stype)){
			trueTree=folder+"TrueTree.tre.";
			for(int i=start; i<end; i++){
				File f = new File(trueTree+i+".pdtree");
				if(f.exists()){
					ReadTree rtree= new ReadTree (trueTree+i+".pdtree");
					double sk=Stat.treeMinBranch(rtree);
					sb.append(i+"	"+sk+"\n");
				}
			}
		}else{
			for(int i=start; i<end; i++){
				File f = new File(trueTree+i+".tre");
				if(f.exists()){
					ReadTree rtree= new ReadTree (trueTree+i+".tre");
					double sk=Stat.treeMinBranch(rtree);
					sb.append(i+"	"+sk+"\n");
				}
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
	
	public void summaryInternalToExternal(String folder, String stype) throws Exception{
		String trueTree=folder+stype+"Tree";
		File compareOut=new File("/lhome/lzhangje/"+stype+".i2e");
		StringBuffer sb=new StringBuffer("");
		if("pda".equals(stype)){
			trueTree=folder+"TrueTree.tre.";
			for(int i=start; i<end; i++){
				File f = new File(trueTree+i+".pdtree");
				if(f.exists()){
					ReadTree rtree= new ReadTree (trueTree+i+".pdtree");
					double sk=Stat.treeInternalExternal(rtree);
					sb.append(i+"	"+sk+"\n");
				}
			}
		}else{
			for(int i=start; i<end; i++){
				File f = new File(trueTree+i+".tre");
				if(f.exists()){
					ReadTree rtree= new ReadTree (trueTree+i+".tre");
					double sk=Stat.treeInternalExternal(rtree);
					sb.append(i+"	"+sk+"\n");
				}
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
	
	
	
	public void summaryIC(String folder, String stype) throws Exception{
		String trueTree=folder+stype+"Tree";
		File compareOut=new File(folder+stype+".ic");
		StringBuffer sb=new StringBuffer("");
		if("pda".equals(stype)){
			trueTree=folder+"TrueTree.tre.";
			for(int i=start; i<end; i++){
				File f = new File(trueTree+i+".pdtree");
				if(f.exists()){
					ReadTree rtree= new ReadTree (trueTree+i+".pdtree");
					double sk=Stat.treeBalance(rtree);
					sb.append(i+"	"+sk+"\n");
				}
			}
		}else{
			for(int i=start; i<end; i++){
				File f = new File(trueTree+i+".tre");
				if(f.exists()){
					ReadTree rtree= new ReadTree (trueTree+i+".tre");
					double sk=Stat.treeBalance(rtree);
					sb.append(i+"	"+sk+"\n");
				}
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
	
	public String getRunTime(String sfin){
		String runtime="no value";
		try{
			File fin=new File(sfin);
			if(fin.exists()){
				FileReader fr=new FileReader(fin);
				BufferedReader bf=new BufferedReader(fr);
				String s=bf.readLine();
				while(s!=null){
					
					if(s.startsWith("Overall execution time")){
						Tokenizer tk=new Tokenizer(s);
						tk.nextToken();
						tk.nextToken();
						tk.nextToken();
						runtime=tk.nextToken();
					}
					
					s=bf.readLine();
				}
			}
		}catch(Exception e){
			e.printStackTrace();
		}
		return runtime;
	}
	
	public void summaryRunTime(String folder, String stype){
		String inforfile=folder+"RAxML_info."+stype;
		
		//File compareOut=new File("/home/zhangje/data/NP/randsample/30_1292.sout");
		File compareOut=new File(folder+stype+".rtime");
		StringBuffer sb=new StringBuffer("");
		
		for(int i=start; i<end; i++){
			String sout=this.getRunTime(inforfile+i);
			if(!"no value".equals(sout)){
				sb.append(i+"	"+sout+"\n");
				System.out.println(i+"	"+sout);
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
	
	
	public double[] jumps (double[] d, double pow){
		double[] dt = new double[d.length+1];
		dt[0]=0;
		double[] jump = new double[d.length];
		for(int i=1; i<d.length; i++){
			dt[i]=Math.pow(d[i-1], -pow);
			System.out.println(d[i-1]);
			System.out.println(dt[i]);
		}
		
		for(int i=0;i<d.length;i++){
			jump[i]=dt[i+1]-dt[i];
		}
		
		return jump;
	}
	
	public void printJump(File fin, File fout){
		double pow=10;
		try{
			FileReader fr=new FileReader(fin);
			BufferedReader bf=new BufferedReader(fr);
			LinkedList l=new LinkedList();
			String s=bf.readLine();
			while(s!=null){
				String[] ss=s.split("	");
				l.add(ss);
				s=bf.readLine();
			}
			bf.close();
			fr.close();
			
			double[] d= new double[l.size()];
			double[] x=new double[l.size()];
			
			for(int i=0;i<l.size();i++){
				String[] ss=(String[])l.get(i);
				x[i]=Double.parseDouble(ss[0]);
				d[i]=Double.parseDouble(ss[1]);
			}
			
			double[] jump=this.jumps(d, pow);
			
			StringBuffer sb = new StringBuffer("");
			
			for(int i=0;i<l.size();i++){
				sb.append(x[i]+"	"+jump[i]+"\n");
			}
			
			
			FileWriter fw=new FileWriter(fout);
			fw.write(sb.toString());
			fw.close();
			
			
		}catch(Exception e){
			e.printStackTrace();
		}
		
		
	}
	
	public double silhouette (double[][] dism, int[] codebook, int[] assignment, int idx_data){
		
		double a= this.avgdissim(dism, idx_data, assignment[idx_data], assignment);
		double b= Double.POSITIVE_INFINITY;
		for(int i=0; i<codebook.length; i++){
			if(assignment[idx_data]!=codebook[i]){
				double tempb=this.avgdissim(dism, idx_data, codebook[i], assignment);
				if (tempb<b){
					b=tempb;
				}
			}
		}
		
		double denominator=a;
		if(b>a){denominator=b;}
		
		return (b-a)/denominator;
	}
	
	public double avgdissim (double[][] dism, int idx_data, int idx_center, int[] assignment){
		double dissim=0.0;
		double count=0;
		for(int i=0;i<assignment.length;i++){
			if(idx_center==assignment[i]){
				dissim+=dism[idx_data][i];
				count+=1;
			}
		}
		return dissim/count;
	}
	
	public double avg_silhouette (double[][] dism, int[] codebook, int[] assignment){
		int numdata=dism[0].length;
		double avg_sil=0.0;
		for(int i=0;i<numdata;i++){
			avg_sil+=this.silhouette(dism, codebook, assignment, i);
		}
		return avg_sil/(double)numdata;
	}
	
	public static void main(String[] args) throws Exception{
		// TODO Auto-generated method stub
		//Analyzer ana =new Analyzer();
		//ana.summaryOutputFiles("/home/zhangje/storage/data/d2554/l1232/PamSample/Pam");
		//ana.summaryRFX("/home/zhangje/storage/data/raxmlXtest/A12/", "WAG");
		//ana.summaryRF("/home/zhangje/storage/data/d2554/l1232/PamSample/", "Pam");
		//ana.summaryRF_NJ("/home/zhangje/storage/data/d2554/l1232/PamSample/", "Pam");
		//ana.summaryRF_PDA("/home/zhangje/storage/data/d2554/l200/PDSample/");
		//ana.summaryRF_NJ_PDA("/home/zhangje/storage/data/d2554/l2000/PDSample/");
		//ana.summaryRF_Par("/home/zhangje/storage/data/d2554/l1000gap/RandSample/", "Rand");
		//System.out.println(ana.getRunTime("/home/zhangje/storage/data/d2554/l200/PhylomapSample/RAxML_info.Phylomap30"));
		//ana.summaryRunTime("/home/zhangje/storage/data/d2554/l200/PDSample/", "pda");
		//ana.summaryIC("/home/zhangje/storage/data/d2554/l1232/PDSample/", "pda");
		
		/*
		ana.summaryMinBranch("/lhome/lzhangje/storage/data/d2554/l1232/PDSample/", "pda");
		ana.summaryMinBranch("/lhome/lzhangje/storage/data/d2554/l1232/PhylomapSample/", "Phylomap");
		ana.summaryMinBranch("/lhome/lzhangje/storage/data/d2554/l1232/SVSample/", "SV");
		ana.summaryMinBranch("/lhome/lzhangje/storage/data/d2554/l1232/CVSample/", "CV");
		ana.summaryMinBranch("/lhome/lzhangje/storage/data/d2554/l1232/RandSample/", "Rand");
		*/
		
		//ana.summaryInternalToExternal("/lhome/lzhangje/storage/data/d2554/l1232/PDSample/", "pda");
		//ana.summaryInternalToExternal("/lhome/lzhangje/storage/data/d2554/l1232/PhylomapSample/", "Phylomap");
		//ana.summaryInternalToExternal("/lhome/lzhangje/storage/data/d2554/l1232/SVSample/", "SV");
		//ana.summaryInternalToExternal("/lhome/lzhangje/storage/data/d2554/l1232/CVSample/", "CV");
		//ana.summaryInternalToExternal("/lhome/lzhangje/storage/data/d2554/l1232/RandSample/", "Rand");
		
		//File fin=new File("/home/zhangje/storage/data/d2554/l1232/Result/SV.err");
		//File fout=new File("/home/zhangje/storage/data/d2554/l1232/Result/SV.err.jump");
		//ana.printJump(fin, fout);
		String truetree = "/lhome/lzhangje/Desktop/RandTree/Par2.tre";
		String infertree = "/lhome/lzhangje/Desktop/RandTree/RAxML_result.ser4x2";
		double dis=TreeDistance.RFdistanceScaled(truetree, infertree);
		System.out.println(dis);
		
	}

}

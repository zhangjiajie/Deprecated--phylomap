package pal.sampling;
import java.util.*;
import java.io.*;

import pal.alignment.*;

public class ColumSample {

	/**
	 * @param args
	 */
	private ReadAlignment aln1;
	private ReadAlignment aln2;
	private String alnnew;
	private double randerr;
	private int laln1;
	private int laln2;
	
	public ColumSample (String aln1, String aln2, String alnnew, double randerr){
		try{
			this.aln1=new ReadAlignment(aln1);
			this.aln2=new ReadAlignment(aln2);
			this.alnnew=alnnew;
			this.randerr=randerr;
			this.laln1=(int)((double)(this.aln1.getSiteCount())*(1-randerr));
			this.laln2=this.aln1.getSiteCount()-this.laln1;
			
		}catch(Exception e){
			e.printStackTrace();
		}
	}
	
	public void toFile (){
		
		Alignment newaln = this.createNewAlignment();
		
		try{
			PrintWriter pw=new PrintWriter(alnnew);
			AlignmentUtils.printInterleaved(newaln, pw);
			pw.close();
		}catch(Exception e){
			e.printStackTrace();
		}
	}
	
	
	private  Alignment createNewAlignment(){
		
		int[] idx=Permutation.next(this.aln1.getSiteCount());
		int[] idx1= new int[this.laln1];
		for(int i=0; i<this.laln1; i++){
			idx1[i]=idx[i];
		}
		
		idx=Permutation.next(this.aln2.getSiteCount());
		int[] idx2= new int[this.laln2];
		for(int i=0; i<this.laln2; i++){
			idx2[i]=idx[i];
		}
		
		ColumSampledAlignment ca1=new ColumSampledAlignment(this.aln1, idx1);
		ColumSampledAlignment ca2=new ColumSampledAlignment(this.aln2, idx2);
		
		Alignment[] alist = {ca1, ca2};
		
		ConcatenatedAlignment alnnew = new ConcatenatedAlignment (alist);
		
		return alnnew;
		
	}
	
	
	public static void main(String[] args) {
		// TODO Auto-generated method stub
		String aln1="/home/zhangje/storage/data/raxmlXtest/A12/WAG1K_TRUE.phy";
		String aln2="/home/zhangje/storage/data/raxmlXtest/A12/ALTER12_TRUE.phy";
		String aln3="/home/zhangje/storage/data/raxmlXtest/A12/MIX_TRUE.phy";
		double randerr=0.2;
		ColumSample cs = new ColumSample(aln1, aln2, aln3, randerr);
		cs.toFile();
	}

}

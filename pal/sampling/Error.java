package pal.sampling;
import pal.alignment.*;
import pal.tree.*;
import pal.distance.*;

import java.io.*;
import java.util.*;

public class Error {

	/**
	 * @param args
	 */
	private double[][] dataMatrix;
	private double[][] sampleDataMatrix;
	private ReadAlignment samplealn;
	private ReadAlignment aln;
	private double avgErr=0;
	private double pd=0;
	
	public Error (String raxmlmatrix, String samplealn, String faln){
		try{
			dataMatrix=DistacneMatrix.readDistacneMatrixRaxml(new File (raxmlmatrix), faln);
			this.samplealn=new ReadAlignment(samplealn);
			this.aln=new ReadAlignment(faln);
			int[] idx = new int[this.samplealn.getSequenceCount()];
			for(int i=0; i<this.samplealn.getSequenceCount(); i++){
				idx[i]=this.aln.whichIdNumber(this.samplealn.getIdentifier(i).getName());
			}
			
			double sumdis=0;
			for(int i=0; i<this.aln.getSequenceCount(); i++){
				double mindis=Double.POSITIVE_INFINITY;
				for(int j=0;j<this.samplealn.getSequenceCount();j++){
					if(mindis<this.dataMatrix[i][idx[j]]){
						mindis=this.dataMatrix[i][idx[j]];
					}
				}
				sumdis+=mindis;
			}
			this.avgErr=sumdis/((double)this.aln.getSequenceCount());
			this.sampleDataMatrix=new double[this.samplealn.getSequenceCount()][this.samplealn.getSequenceCount()];
			for(int i=0; i<this.samplealn.getSequenceCount(); i++){
				for (int j=0; j<this.samplealn.getSequenceCount(); j++){
					this.sampleDataMatrix[i][j]=this.dataMatrix[idx[i]][idx[j]];
				}
			}
		}catch(Exception e){
			e.printStackTrace();
		}
	}
	
	public Error (double[][] dataMatrix, String samplealn, ReadAlignment aln){
		try{
			this.dataMatrix=dataMatrix;
			this.samplealn=new ReadAlignment(samplealn);
			this.aln=aln;
			int[] idx = new int[this.samplealn.getSequenceCount()];
			for(int i=0; i<this.samplealn.getSequenceCount(); i++){
				idx[i]=this.aln.whichIdNumber(this.samplealn.getIdentifier(i).getName());
			}
			this.sampleDataMatrix=new double[this.samplealn.getSequenceCount()][this.samplealn.getSequenceCount()];
			double sumdis=0;
			for(int i=0; i<this.aln.getSequenceCount(); i++){
				double mindis=Double.POSITIVE_INFINITY;
				for(int j=0;j<this.samplealn.getSequenceCount();j++){
					if(this.dataMatrix[i][idx[j]]<mindis){
						mindis=this.dataMatrix[i][idx[j]];
					}
				}
				sumdis+=mindis;
			}
			this.avgErr=sumdis/((double)this.aln.getSequenceCount());
			
			for(int i=0; i<this.samplealn.getSequenceCount(); i++){
				for (int j=0; j<this.samplealn.getSequenceCount(); j++){
					this.sampleDataMatrix[i][j]=this.dataMatrix[idx[i]][idx[j]];
				}
			}
			
		}catch(Exception e){
			e.printStackTrace();
		}
	}
	
	public Error (String raxmlmatrix, String samplealn, ReadAlignment aln){
		try{
			dataMatrix=DistacneMatrix.readDistacneMatrixRaxml(new File (raxmlmatrix), aln);
			this.samplealn=new ReadAlignment(samplealn);
			this.aln=aln;
			int[] idx = new int[this.samplealn.getSequenceCount()];
			for(int i=0; i<this.samplealn.getSequenceCount(); i++){
				idx[i]=this.aln.whichIdNumber(this.samplealn.getIdentifier(i).getName());
			}
			
			double sumdis=0;
			for(int i=0; i<this.aln.getSequenceCount(); i++){
				double mindis=Double.POSITIVE_INFINITY;
				for(int j=0;j<this.samplealn.getSequenceCount();j++){
					if(this.dataMatrix[i][idx[j]]<mindis){
						mindis=this.dataMatrix[i][idx[j]];
					}
				}
				sumdis+=mindis;
			}
			this.avgErr=sumdis/((double)this.aln.getSequenceCount());
			this.sampleDataMatrix=new double[this.samplealn.getSequenceCount()][this.samplealn.getSequenceCount()];
			for(int i=0; i<this.samplealn.getSequenceCount(); i++){
				for (int j=0; j<this.samplealn.getSequenceCount(); j++){
					this.sampleDataMatrix[i][j]=this.dataMatrix[idx[i]][idx[j]];
				}
			}
		}catch(Exception e){
			e.printStackTrace();
		}
	}
	
	public Error(){
		
	}
	
	public double getAverageError(){
		return this.avgErr;
	}
	
	private void pdcal(Node root){
		for(int i=0; i<root.getChildCount();i++){
			this.pdcal(root.getChild(i));
		}
		this.pd+=root.getBranchLength();
	}
	
	public double getPDvalue(SimpleTree st){
		this.pdcal(st.getRoot());
		return this.pd;
	}
	
	public DistanceMatrix getSampleDistanceMatrix(){
		return new DistanceMatrix(this.sampleDataMatrix, this.samplealn.getIdGroup());
	}
	
	
	public static void main(String[] args) throws Exception {
		// TODO Auto-generated method stub
		String faln="/home/zhangje/storage/data/d2554/l200/alignment_sim.phy";
		ReadAlignment aln = new ReadAlignment(faln);
		double[][] datamatrix = DistacneMatrix.readDistacneMatrixRaxml(new File ("/home/zhangje/storage/data/d2554/l200/raxmlmatrix.mx"), faln);
		DistanceMatrix dm= new DistanceMatrix(datamatrix, aln.getIdGroup());
		NeighborJoiningTree njt = new NeighborJoiningTree (dm);
		ExtractSubtree.writeSubTree(new File("/home/zhangje/storage/data/d2554/l200/NJ.tre"), njt);
	}

}

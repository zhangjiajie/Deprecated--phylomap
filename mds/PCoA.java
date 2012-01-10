package mds;
import java.util.*;
import java.io.*;

import javastat.multivariate.*;

import Jama.*;

public class PCoA implements Serializable{

	/**
	 * @param args
	 */
	public double[] Varience;
	public double[][] distanceMatrix;
	
	//sorted eignvectors and value
	public double[][] eigenVector;
	public double[] eigenValue;
	
	public double[][] getDistanceMatrix(){
		return this.distanceMatrix;
	}
	
	public List readAlignedSeq(File fin){
		LinkedList l=new LinkedList();
		
		try{
			
		}catch(Exception e){
			e.printStackTrace();
		}
		return l;
	}
	
	public double dotProduct(double[] x, double[] y){
		double acc=0;
		for(int i=0;i<x.length;i++){
			
			acc=acc+x[i]*y[i];
			
		}
		
		return acc;
	}
	
	public double[] getColum(double[][] dm, int k){
		double[] colum=new double[dm.length];
		
		for(int i=0;i<dm.length;i++){
			colum[i]=dm[i][k];
		}
		
		return colum;
	}
	
	public double[][] transpose(double[][] M){
		int nrow=M.length;
		int ncol=M[0].length;
		double[][] Mt=new double[ncol][nrow];
		for(int i=0;i<ncol;i++){
			Mt[i]=this.getColum(M, i);
		}
		return Mt;
	}
	
	public void toFile(File fout, double[][] coords){
		try{
			StringBuffer sb=new StringBuffer("");
			int m=coords.length;
			int n=coords[0].length;
			for(int i=0;i<m;i++){
				for(int j=0;j<n;j++){
					sb.append(coords[i][j]+"	");
				}
				sb.append("\n");
			}
			
			FileWriter fw=new FileWriter(fout);
			fw.write(sb.toString());
			fw.close();
		}catch (Exception e){
			e.printStackTrace();
		}
	}
	
	public double sumSquare(double[] row){
		double mean=0;
		int n=row.length;
		for(int i=0;i<n;i++){
			mean=mean+row[i]*row[i];
		}
		
		return mean;
	}
	
	public double rowMean(double[] row){
		double mean=0;
		int n=row.length;
		for(int i=0;i<n;i++){
			mean=mean+row[i];
		}
		mean=mean/((double)n);
		return mean;
	}
	
	public double matrixMean(double[][] M){
		
		double mean=0;
		int m=M.length;
		int n=M[0].length;
		for(int i=0;i<m;i++){
			for(int j=0;j<n;j++){
				mean=mean+M[i][j];
			}
		}
		mean=mean/((double)(m*n));
		return mean;
	}
	
	public double[][] centrolizeMatrix(double[][] M){
		int n=M.length;
		double[][] E=new double[n][n];
		double[][] F=new double[n][n];
		for(int i=0;i<n;i++){
			for(int j=0;j<n;j++){
				E[i][j]=M[i][j]*M[i][j]*(-1.0/2.0);
			}
		}
		double e_mean=this.matrixMean(E);
		double[] ej_=new double[n];
		double[] e_k=new double[n];
		for(int i=0;i<n;i++){
			ej_[i]=this.rowMean(E[i]);
			double[] col=this.getColum(E, i);
			e_k[i]=this.rowMean(col);
		}
		for(int i=0;i<n;i++){
			for(int j=0;j<n;j++){
				F[i][j]=E[i][j]-ej_[i]-e_k[j]+e_mean;
			}
		}
		return F;
	}
	
	// read n x n matrix from phylip calculation
	public double[][] PcoordPhylip(double[][] disM){
		double[][] disMatrix=disM;
		this.distanceMatrix=disMatrix;
		double[][] f=this.centrolizeMatrix(disMatrix);
		Matrix F=new Matrix(f);
		EigenvalueDecomposition ed=new EigenvalueDecomposition(F);
		double[] eigenValue=ed.getRealEigenvalues();
		Matrix eigenM=ed.getV();
		double[][] eM=eigenM.transpose().getArray();
		int n=f.length;
		// Ordering
		boolean isChange=true;
		while(isChange){
			isChange=false;
			for(int i=0;i<n-1;i++){
				if(eigenValue[i+1]>eigenValue[i]){
					double temp=eigenValue[i+1];
					eigenValue[i+1]=eigenValue[i];
					eigenValue[i]=temp;
					double[] atemp=eM[i+1];
					eM[i+1]=eM[i];
					eM[i]=atemp;
					isChange=true;
				}
			}
		}
		
		this.eigenVector=eM;
		this.eigenValue=eigenValue;
		
		double [] firstPC=eM[0];
		double [] secondPC=eM[1];
		double[][] pcaResult=new double[3][n];
		pcaResult[0]=firstPC;
		pcaResult[1]=secondPC;
		pcaResult[2]=eigenValue;
		pcaResult=this.normalize(pcaResult);
		return pcaResult;
	}
	
	//write the accumulative eignvalues to file for plotting
	public double accumulateEigenvaluesTofile(File fout){
		StringBuffer sb=new StringBuffer("");
		double feigen=0;
		if(this.eigenValue!=null){
			int len=this.eigenValue.length;
			double sum=0;
			for(int i=0;i<len;i++){
				if(this.eigenValue[i]>0){
					sum=sum+this.eigenValue[i];
				}
			}
			double accum=0;
			feigen=this.eigenValue[0]/sum;
			for(int i=0;i<len;i++){
				if(this.eigenValue[i]>0){
					accum=accum+this.eigenValue[i];
					double ration=accum/sum;
					sb.append((i+1)+"	"+ration+"\n");
				}
			}
			
			try{
				FileWriter fw=new FileWriter(fout);
				fw.write(sb.toString());
				fw.close();
			}catch(Exception e){
				e.printStackTrace();
			}
			
		}
		
		return feigen;
	}
	
	//return the first m eigen vectors whose accumulative eigenvalues exceed p percent of total
	public double[][] getFirstMdimension(double p){
		
		if(this.eigenValue!=null){
			int len=this.eigenValue.length;
			int n=this.eigenVector[0].length;
			double sum=0;
			for(int i=0;i<len;i++){
				if(this.eigenValue[i]>0){
					sum=sum+this.eigenValue[i];
				}
			}
			double accum=0;
			int cnt=0;
			for(int i=0;i<len;i++){
				if(this.eigenValue[i]>0){
					accum=accum+this.eigenValue[i];
					double ration=accum/sum;
					if(ration<p){
						cnt++;
					}
				}
			}
			double[][] eigenVm=new double[cnt][n];
			for(int i=0;i<cnt;i++){
				eigenVm[i]=this.eigenVector[i];
			}
			return eigenVm;
		}else{
			return null;
		}
	}
	
	public double[][] normalize(double[][] coords){
		
		double k0=this.sumSquare(coords[0]);
		double k1=this.sumSquare(coords[1]);
		k0=coords[2][0]/k0;
		k1=coords[2][1]/k1;
		k0=Math.sqrt(k0);
		k1=Math.sqrt(k1);
		for(int i=0;i<coords[0].length;i++){
			coords[0][i]=k0*coords[0][i];
			coords[1][i]=k1*coords[1][i];
		}
		return coords;
	}
	
	public double[][] normalizeHD(double[][] coords){
		int n=coords[0].length;
		
		for(int k=0;k<n;k++){
			double k0=this.sumSquare(coords[k]);
			k0=coords[n][k]/k0;
			k0=Math.sqrt(k0);
			for(int i=0;i<coords[0].length;i++){
				coords[k][i]=k0*coords[k][i];
				
			}
		}
		
		return coords;
	}
	
	public double[][] normalize2(double[][] coords){
		
		double k0=Math.sqrt(coords[2][0]);
		double k1=Math.sqrt(coords[2][1]);
		for(int i=0;i<coords[0].length;i++){
			coords[0][i]=coords[0][i]/k0;
			coords[1][i]=coords[1][i]/k1;
		}
		return coords;
	}
	
	public static void main(String[] args) {
		// TODO Auto-generated method stub

	}

}

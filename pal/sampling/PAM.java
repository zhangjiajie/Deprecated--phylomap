package pal.sampling;
import java.io.*;
import java.util.*;

public class PAM {

	/**
	 * @param args
	 */
	
	double[][] distanceMatrix;
	int numClusters=1;
	int[] codebook;
	int[] finalAssignment;	
	int[] unseleted;
	int numData=0;
	int maxIters=numData*1000;
	double finalcost=0.0;
	
	public PAM(String FdistanceMatrix, int numClusters){
		
	}
	
	public PAM(double[][] distanceMatrix, int numClusters){
		this.distanceMatrix=distanceMatrix;
		this.numData=this.distanceMatrix[0].length;
		this.numClusters=numClusters;
		this.codebook=new int[this.numClusters];
		this.unseleted=new int[this.numData-this.numClusters];
		this.maxIters=numData*1000;
		this.learn();
	}
	
	private int[] assign(){
		int[] assignment = new int[numData];
		for(int i=0;i<numData; i++){
			double shortdis=distanceMatrix[i][codebook[0]];
			int idx=0;
			for(int j=1;j<numClusters;j++){
				if(distanceMatrix[i][codebook[j]]< shortdis){
					shortdis=distanceMatrix[i][codebook[j]];
					idx=j;
				}
			}
			assignment[i]=codebook[idx];
			//System.out.println(assignment[i]);
		}
		return assignment;
	}
	
	private double cost(int[] assignment){
		
		double cost=0.0;
		for(int i=0;i<numData;i++){
			cost+=this.distanceMatrix[i][assignment[i]];
		}
		return cost;
	}
	
	
	private void build(){
		int[] randidx=Permutation.next(this.numData);
		for(int i=0;i<numClusters;i++){
			codebook[i]=randidx[i];
		}
		for(int i=0; i<(numData-numClusters);i++){
			unseleted[i]=randidx[i+numClusters];
		}
	}
	
	private boolean swap(int[] assignment){
		boolean changed=false;
		for(int i=0;i<numClusters;i++){
			double mincost= this.cost(assignment);
			for(int j=0; j< (numData-numClusters);j++){
				int oldidx=codebook[i];
				int newidx=unseleted[j];
				codebook[i]=newidx;
				unseleted[j]=oldidx;
				int[] newassignment=this.assign();
				double newcost=this.cost(newassignment);
				if(newcost>=mincost){
					codebook[i]=oldidx;
					unseleted[j]=newidx;
				}else{
					mincost=newcost;
					changed=true;
				}
			}
			
		}
		return changed;
	}
	
	public void learn(){
		build();
		boolean changed=true;
		int count=0;
		int[] assignment=null;
		while(changed && count<maxIters){
			changed=false;
			count++;
			assignment=this.assign();
			changed=swap(assignment);
			//System.out.println("Iters:"+count);
		}
		this.finalAssignment=assignment;
		//this.finalcost=this.cost(this.finalAssignment);
	}
	
	public int[] getCodebook(){
		return this.codebook;
	}
	
	public int[] getAssignment(){
		return this.finalAssignment;
	}
	
	
	
	public static void main(String[] args) throws Exception{
		// TODO Auto-generated method stub
			double[][] dataMatrix=DistacneMatrix.readDistacneMatrixRaxml(new File("/home/zhangje/NewData/raxmlmatrix.mx"), "/home/zhangje/NewData/alignment_sim.phy");
			double[] avg_s=new double[120];
			for(int i=80;i<200;i++){
				PAM pam = new PAM(dataMatrix, i);
				int[] codebook = pam.getCodebook();
				int[] assignment = pam.getAssignment();
				Analyzer ana = new Analyzer ();
				avg_s[i-80]=ana.avg_silhouette(dataMatrix, codebook, assignment);
				System.out.println(i+": "+avg_s[i-80]);
			}
		
	}

}

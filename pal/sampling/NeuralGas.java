package pal.sampling;

import java.util.*;
import java.io.*;

public class NeuralGas {

	/**
	 * @param args
	 */
	int num_of_cluster=10;
	int tmax=10;
	double epsilon_start=1;
	double epsilon_end=0.001;
	double ramda_start=num_of_cluster;
	double ramda_end=0.001;
	double[][] data; // input data matrix
	double[][] centers; // the values of data centers
	int N=0; //number of data points
	int d=0; //dim
	int[] codebook; // the index list of data center
	int[] cluster_idx; // each data point belongs to which cluster center
	double error=0;
	
	public NeuralGas(){
		
	}
	
	public NeuralGas(double[][] data, int num_of_cluster, int tmax){
		this.data=data;
		this.num_of_cluster=num_of_cluster;
		this.ramda_start=num_of_cluster;
		this.tmax=tmax;
		N=data.length;
		d=data[0].length;
		codebook=new int[num_of_cluster];
		centers=new double[num_of_cluster][d];
		cluster_idx=new int[N];
	}
	
	public static double norm2vector(double[] v1, double[] v2){
		double dis=0;
		for(int i=0;i<v1.length;i++){
			dis=dis+(v1[i]-v2[i])*(v1[i]-v2[i]);
		}
		return dis;
	}
	
	//v1+v2
	public static double[] vector_add(double[] v1, double[] v2){
		
		double[] v=new double[v1.length];
		for(int i=0;i<v1.length;i++){
			v[i]=v1[i]+v2[i];
		}
		return v;
	}
	
	//v1-v2
	public static double[] vector_minus(double[] v1, double[] v2){
		
		double[] v=new double[v1.length];
		for(int i=0;i<v1.length;i++){
			v[i]=v1[i]-v2[i];
		}
		return v;
	}
	
	//f*v1
	public static double[] vector_factor_mutilple(double[] v1, double f){
		
		double[] v=new double[v1.length];
		for(int i=0;i<v1.length;i++){
			v[i]=v1[i]*f;
		}
		return v;
	}
	
	
	public void learn(){
		//Random init num_of_cluster
		Date data1=new Date();
		Random rand=new Random(data1.getTime());
		HashSet hs=new HashSet();
		for(int i=0;i<this.num_of_cluster;i++){
			int tmp=rand.nextInt(N);
			Integer I=new Integer(tmp);
			while(hs.contains(I)){
				tmp=rand.nextInt(N);
				I=new Integer(tmp);
			}
			hs.add(I);
			codebook[i]=tmp;
			centers[i]=data[tmp];
		}
		
		double epsilon=0;
		double ramda=0;
		double[] dists=new double[num_of_cluster];
		double[] rank=new double[num_of_cluster];
		HashSet hsstop = new HashSet ();
		
		//Start pattern by pattern learning 
		for(int i=0;i<tmax;i++){
			System.out.println("Epoc: "+i);
			//Cooling down the learning rate 
			epsilon=epsilon_start*(Math.pow((epsilon_end/epsilon_start), ((double)i/(double)tmax)));
			
			//Narrow down the neighborhood 
			ramda=ramda_start*(Math.pow((ramda_end/ramda_start), ((double)i/(double)tmax)));
			
			//Learn each point but in random order
			int[] idx=Permutation.next(N);
			for(int j=0;j<N;j++){
				// get the distance from current point to all cluster centers 
				double[] currentPoint=data[idx[j]];
				for(int k=0;k<num_of_cluster;k++){
					dists[k]=norm2vector(currentPoint, centers[k]);
				}
				
				// get the rank rank[i]=number of center[j] with norm(center[j],currentPoint)<norm(center[i],currentPoint)
				int[] IX=this.order(dists);
				for(int k=0;k<num_of_cluster;k++){
					rank[IX[k]]=k;
				}
				
				// update all cluster centers according to the rank
				for(int k=0;k<num_of_cluster;k++){
					double factor=epsilon*Math.exp((-rank[k]/ramda));
					double[] v=vector_minus(currentPoint, centers[k]);
					v=vector_factor_mutilple(v, factor);
					centers[k]=vector_add(centers[k],v);
				}
			}
			
			boolean hasChanged = false;
			HashSet hsstop2=new HashSet();
			for(int m=0;m<num_of_cluster;m++){
				double mindis=Double.POSITIVE_INFINITY;
				int idxm=-1;
				for(int j=0;j<N;j++){
					double[] currentPoint=data[j];
					double dis=norm2vector(currentPoint, centers[m]);
					if(dis<mindis){
						mindis=dis;
						idxm=j;
					}
				}
				//centers[m]=data[idxm];
				if(hsstop.contains(idxm+"")){
					hsstop2.add(idxm+"");
				}else{
					hsstop2.add(idxm+"");
					hasChanged = true;
				}
			}
			hsstop=hsstop2;
			
			if(!hasChanged){
				break;
			}
			
			
		}
		
		//Find the closest data point to the cluster center codebook, modified to avoid duplication
		for(int i=0;i<num_of_cluster;i++){
			double mindis=Double.POSITIVE_INFINITY;
			int idx=-1;
			for(int j=0;j<N;j++){
				double[] currentPoint=data[j];
				double dis=norm2vector(currentPoint, centers[i]);
				if(dis<mindis){
					boolean incb=false;
					for(int k=0; k<i; k++){
						if(codebook[k]==j){
							incb=true;
						}
					}
					if(!incb){
						mindis=dis;
						idx=j;
					}
				}
			}
			codebook[i]=idx;
		}
		
		//Assign each data point to its center
		for(int i=0;i<N;i++){
			double[] currentPoint=data[i];
			double mindis=Double.POSITIVE_INFINITY;
			int idx=-1;
			for(int j=0;j<num_of_cluster;j++){
				double dis=norm2vector(currentPoint, data[codebook[j]]);
				if(dis<mindis){
					mindis=dis;
					idx=j;
				}
			}
			cluster_idx[i]=codebook[idx];
		}
		
		
		//Get the quantization error
		for(int i=0;i<N;i++){
			error=error+norm2vector(data[i], data[cluster_idx[i]]);
		}
		error=error/(double)N;
		
		
	}
	
	// ascending order return the index
	public int[] order(double[] A){
		int[] Index=new int[A.length];
		for(int i=0;i<A.length;i++){
			Index[i]=i;
		}
		boolean isChange=true;
		while(isChange){
			isChange=false;
			for(int i=0;i<A.length-1;i++){
				if(A[i+1]<A[i]){
					double temp=A[i+1];
					A[i+1]=A[i];
					A[i]=temp;
					int atemp=Index[i+1];
					Index[i+1]=Index[i];
					Index[i]=atemp;
					isChange=true;
				}
			}
			
			
		}
		
		return Index;
	}
	
	public int[] getCodebookvector(){
		return this.codebook;
	}
	
	public int[] getClusterindex(){
		return this.cluster_idx;
	}
	
	public double getError(){
		return this.error;
	}
	
	
	public static void main(String[] args) {
		// TODO Auto-generated method stub
		double[] A={5,4,3,2,1,0};
		NeuralGas ng= new NeuralGas();
		int[] IX=ng.order(A);
		for(int i=0;i<IX.length;i++){
			//System.out.println(IX[i]);
		}
		
		File f=new File("/home/zhangjiajie/Desktop/data");
		try{
			ArrayList al=new ArrayList();
			FileReader fr=new FileReader(f);
			BufferedReader bf=new BufferedReader(fr);
			String s=bf.readLine();
			while(s!=null){
				al.add(s.trim());
				s=bf.readLine();
			}
			double[][] D=new double[al.size()][2];
			bf.close();
			fr.close();
			
			for(int i=0;i<al.size();i++){
				String[] ss=((String)al.get(i)).split("	");
				D[i][0]=Double.parseDouble(ss[0]);
				D[i][1]=Double.parseDouble(ss[1]);
			}
			
			ng=new NeuralGas(D,5,10);
			ng.learn();
			int[] codebook=ng.getCodebookvector();
			StringBuffer sb=new StringBuffer("");
			for(int i=0;i<codebook.length;i++){
				sb.append(D[codebook[i]][0]+"	"+D[codebook[i]][1]+"\n");
				
			}
			File fout=new File("/home/zhangjiajie/Desktop/center.dat");
			FileWriter fw3=new FileWriter(fout);
			fw3.write(sb.toString());
			fw3.close();
			
			
			for(int i=1;i<100;i++){
				ng=new NeuralGas(D,i,10);
				ng.learn();
				System.out.println(ng.getError());
			}
			
		}catch(Exception e){
			e.printStackTrace();
		}
		
	}

}

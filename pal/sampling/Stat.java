package pal.sampling;

import static java.lang.Math.pow;
import pal.tree.*;

import java.util.*;


public class Stat {

	/**
	 * @param args
	 */
	
    /**
     * Calculates the sum of the data.
     * @param data the input data (two dimensional array).
     * @return the sum of the data.
     */

    public static double sum(double[] ...data)
    {
        double s = 0;
        for (int i = 0; i < data.length; i++)
        {
            for (int j = 0; j < data[i].length; j++)
            {
                s += data[i][j];
            }
        }

        return s;
    }

    public static double sum2(double[] ...data)
    {
        double s = 0;
        for (int i = 0; i < data.length; i++)
        {
            for (int j = 0; j < data[i].length; j++)
            {
                s += (data[i][j]*data[i][j]);
            }
        }

        return s;
    }
    /**
     * Calculates the mean of the input data.
     * @param data the input data.
     * @return the mean of the input data.
     * @exception IllegalArgumentException the length of the input data should
     *                                     not be 0.
     */

    public static double mean(double[] data)
    {
        if (data.length == 0)
        {
            throw new IllegalArgumentException(
                    "The length of the input data should not be 0.");
        }

        return sum(data) / data.length;
    }
	
    public static double mean2(double[] data)
    {
        if (data.length == 0)
        {
            throw new IllegalArgumentException(
                    "The length of the input data should not be 0.");
        }
       
        return sum2(data) / data.length;
    }
    
    /**
     * Calculates the variance of the input data.
     * @param data the input data.
     * @return the variance of the input data.
     * @exception IllegalArgumentException the length of the input data should
     *                                     not be 0.
     */

    public static double variance(double[] data)
    {
        if (data.length == 0)
        {
            throw new IllegalArgumentException(
                    "The length of the input data should not be 0.");
        }
        double var=0;
        
        var=mean2(data)-pow (mean(data),2);
        
        return var; 
    }
    
    public static double sd(double[] data){
    	
    	return Math.sqrt(variance (data));
    }
    
    public static double skewness(double[] data){
    	double sk=0;
    	
    	double mean=mean(data);
    	double sd=sd(data);
    	
    	 for (int i = 0; i < data.length; i++)
         {
    		 sk+=pow((data[i]-mean)/sd,3);  
         }
    	
    	return sk/data.length;
    	
    }
    
    public static double min (double[] data){
    	double min = 0;
    	/*
    	for(int i=0; i<data.length; i++){
    		if(data[i] < min){
    			min =data[i];
    		}
    	}
    	*/
    	//double sd = sd (data);
    	double mean = mean (data);
    	//sd = mean - 0.5 * sd ;
    	//System.out.println(sd);
    	min =0 ;
    	Arrays.sort(data);
    	
    	
    	int L = (data.length / 10); 
    	for(int i=0; i< L ; i++){
    		if(data[i] < mean){
    			min += data[i];
    		}else{
    			L--;
    		}
    	}
    	
    	/*
    	int cnt =0;
    	for(int i=0; i<data.length; i++){
    		if(data[i] < sd){
    			min += data[i];
    			cnt++;
    		}
    	}
    	*/
    	return min / (double) L;
    }
    
    public static double treeSkewness (Tree t){
    	
    	LinkedList<Double> blist=new LinkedList<Double>();
    	for(int i=0 ;i<t.getInternalNodeCount();i++){
    		Node n=t.getInternalNode(i);
    		if(!n.isRoot()){
    			blist.add(n.getBranchLength());
    		}
    	}
    	
    	for(int i=0;i<t.getExternalNodeCount();i++){
    		Node n=t.getExternalNode(i);
    		if(n.getBranchLength() > 0){
    			blist.add(n.getBranchLength());
    		}
    	}
    	double[] data = new double[blist.size()];
    	for(int i=0; i<blist.size(); i++){
    		data[i]=blist.get(i);
    	}
    	
    	return skewness(data);
    }
    
    public static double treeVariance (Tree t){
    	LinkedList<Double> blist=new LinkedList<Double>();
    	for(int i=0 ;i<t.getInternalNodeCount();i++){
    		Node n=t.getInternalNode(i);
    		if(!n.isRoot()){
    			blist.add(n.getBranchLength());
    		}
    	}
    	
    	for(int i=0;i<t.getExternalNodeCount();i++){
    		Node n=t.getExternalNode(i);
    		if(n.getBranchLength() > 0){
    			blist.add(n.getBranchLength());
    		}
    	}
    	double[] data = new double[blist.size()];
    	for(int i=0; i<blist.size(); i++){
    		data[i]=blist.get(i);
    	}
    	
    	return variance(data);
    }
    
    public static double treeMinBranch (Tree t){
    	LinkedList<Double> blist=new LinkedList<Double>();
    	for(int i=0 ;i<t.getInternalNodeCount();i++){
    		Node n=t.getInternalNode(i);
    		if(!n.isRoot()){
    			if(n.getBranchLength() > 0){
    				blist.add(n.getBranchLength());
    			}
    		}
    	}
    	
    	for(int i=0;i<t.getExternalNodeCount();i++){
    		Node n=t.getExternalNode(i);
    		if(n.getBranchLength() > 0){
    			blist.add(n.getBranchLength());
    		}
    	}
    	double[] data = new double[blist.size()];
    	for(int i=0; i<blist.size(); i++){
    		data[i]=blist.get(i);
    	}
    	
    	return min(data);
    }
    
    
    public static double treeInternalExternal (Tree t){
    	LinkedList<Double> blist=new LinkedList<Double>();
    	LinkedList<Double> alist=new LinkedList<Double>();
    	for(int i=0 ;i<t.getInternalNodeCount();i++){
    		Node n=t.getInternalNode(i);
    		if(!n.isRoot()){
    			//if(n.getBranchLength() > 0){
    				alist.add(n.getBranchLength());
    			//}
    		}
    	}
    	
    	for(int i=0;i<t.getExternalNodeCount();i++){
    		Node n=t.getExternalNode(i);
    		//if(n.getBranchLength() > 0){
    			blist.add(n.getBranchLength());
    		//}
    	}
    	double[] dataa = new double[alist.size()];
    	double[] datab = new double[blist.size()];
    	for(int i=0; i<blist.size(); i++){
    		datab[i]=blist.get(i);
    	}
    	
    	for(int i=0; i<alist.size(); i++){
    		dataa[i]=alist.get(i);
    	}
    	return mean(dataa) / mean(datab);
    }
    
    
    public static double treeBalance (Tree t){
    	double IC=0.0;
    	
    	for(int i=0; i<t.getInternalNodeCount(); i++){
    		//System.out.println(i);
    		SimpleNode node=(SimpleNode)t.getInternalNode(i);
    		//System.out.println(node.getNumber());
    		
    		Stat s1 = new Stat();
    		s1.getLeafCnt(node.getChild(0));
    		
    		Stat s2= new Stat();
    		s2.getLeafCnt(node.getChild(1));
    		
    		IC+=Math.abs(s1.leafcnt-s2.leafcnt);
    		//System.out.println("tl:"+s1.leafcnt);
    		
    		//System.out.println("tr:"+s2.leafcnt);
    	}
    	
    	
    	
    	double n=t.getExternalNodeCount();
    	
    	return IC/((n-1)*(n-2)/2.0);
    }
    
    public int leafcnt=0;
    
    public void getLeafCnt(Node n){
    	
    	if(n.isLeaf()){
    		leafcnt++;
    		return;
    	}else{
    		for(int i=0; i<n.getChildCount(); i++){
    			this.getLeafCnt(n.getChild(i));
    		}
    	}
    }
	
	public static void main(String[] args) throws Exception{
		// TODO Auto-generated method stub
		ReadTree rtree= new ReadTree ("/home/zhangje/storage/data/d2554/l1232/TrueTree.tre");
		System.out.println(Stat.treeBalance(new SimpleTree(rtree.getRoot())));
	}

}

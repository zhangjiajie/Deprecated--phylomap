/**
 * PNP: PCA in assosionation with Neural gas sampling and Phlogentic tree representation
 */
package phylomap;
import java.util.*;
import java.io.*;

/**
 * @author zhangjiajie
 *
 */
public class postOrder {

	/**
	 * @param args
	 */
	
	// nodename to 3 connected nodes distances 
	HashMap hm_innernodes_dis=new HashMap();
	
	
	HashMap hm_mds_innernodes_dis=new HashMap();
	
	//In use 
	//nodeA:nodeB distance, store distances between 2 nodes
	HashMap hm_innernodes_dis2=new HashMap();
	
	//In use 
	//distance matrix for shannon mapping 
	HashMap hm_mds_innernodes_dis2=new HashMap();
	
	HashMap hm_node_connectednodename=new HashMap();
	
	//All nodes inner and leaf nodes list
	ArrayList nodes_list=new ArrayList();
	
	HashMap hm_controlpoint=new HashMap();
	pnp PNP=new pnp();
	Visual vis=new Visual();
	
	StringBuffer sb=new StringBuffer("");
	
	double mf=0.1;
	
	double scaleFactor=1;
	
	double error=0.0001;
	
	int counter=0;
	
	double smallestdis=Double.POSITIVE_INFINITY;
	
	//In use
	public ArrayList getNodeList(){
		return this.nodes_list;
	}
	
	//In use 
	public void setScaleFactor(double scale){
		this.scaleFactor=scale;
	}
	
	
	public double dis(double[] c1, double[] c2){
		double dis=0;
		for(int i=0;i<c1.length;i++){
			dis+=((c1[i]-c2[i])*(c1[i]-c2[i]));
		}
		return Math.sqrt(dis);
	}
	
	//In use 
	//Find the smallest distance in the tree that > 0
	public void postTreefindthesmallestDistance(node n){
		if(n.getLeft()!=null && n.getRight()!=null){
			postTreefindthesmallestDistance(n.getLeft());
			postTreefindthesmallestDistance(n.getRight());
			if(n.getDistance()<this.smallestdis && n.getDistance()>0){
				this.smallestdis=n.getDistance();
			}
			
		}else{
			if(n.getDistance()<this.smallestdis && n.getDistance()>0){
				this.smallestdis=n.getDistance();
			}
			
			return;
		}
	}
	
	
	//In use 
	//Correct the minus distance 
	public void postTreeCorrectMinusDis(node n){
		if(n.getLeft()!=null && n.getRight()!=null){
			postTreeCorrectMinusDis(n.getLeft());
			postTreeCorrectMinusDis(n.getRight());
			if(n.getDistance()<=0){
				n.setDistance(this.smallestdis);
			}
			
		}else{
			if(n.getDistance()<=0){
				n.setDistance(this.smallestdis);
			}
			
			return;
		}
	}
	
	//In use 
	//Correct the minus distance 
	public void correctMinusDisInTree(node root){
		this.postTreefindthesmallestDistance(root.getLeft());
		this.postTreeCorrectMinusDis(root.getLeft());
	}
	
	//In use 
	//put all nodes name into nodes_list
	public void postTree(node n){
		if(n.getLeft()!=null && n.getRight()!=null){
			postTree(n.getLeft());
			postTree(n.getRight());
			counter++;
			this.nodes_list.add(n.getId());
		}else{
			counter++;
			this.nodes_list.add(n.getId());
			return;
		}
	}
	
	//In use
	//Fill in hm_innernodes_dis: nodename to 3 connected nodes distances  
	public void postTree_cal_tree_dis(node n){
		if(n.getLeft()!=null && n.getRight()!=null){
			postTree_cal_tree_dis(n.getLeft());
			postTree_cal_tree_dis(n.getRight());
			double[] X=new double[3];
			X[0]=n.getDistance();
			X[1]=n.getLeft().getDistance();
			X[2]=n.getRight().getDistance();			
			hm_innernodes_dis.put(n.getId().trim(), X);
		}else{
			return;
		}
	}
	
	//In use
	//Fill in hm_innernodes_dis2:  nodeA:nodeB distance
	public void postTree_cal_tree_dis2(node n){
		if(n.getLeft()!=null && n.getRight()!=null){
			postTree_cal_tree_dis2(n.getLeft());
			postTree_cal_tree_dis2(n.getRight());
			double[] X=new double[3];
			X[0]=n.getDistance()*scaleFactor;
			X[1]=n.getLeft().getDistance()*scaleFactor;
			X[2]=n.getRight().getDistance()*scaleFactor;
			String namec=n.getId().trim();
			String namel=n.getLeft().getId().trim();
			String namer=n.getRight().getId().trim();
			hm_innernodes_dis2.put(namec+":"+namel, new Double(X[1]));
			hm_innernodes_dis2.put(namel+":"+namec, new Double(X[1]));
			hm_innernodes_dis2.put(namec+":"+namer, new Double(X[2]));
			hm_innernodes_dis2.put(namer+":"+namec, new Double(X[2]));
		}else{
			return;
		}
	}
	
	
	
	
	public void postTree_cal_mds_dis(node n, HashMap hm){
		if(n.getLeft()!=null && n.getRight()!=null){
			postTree_cal_mds_dis(n.getLeft(),hm);
			postTree_cal_mds_dis(n.getRight(),hm);
			//System.out.println(n.getId()+":"+n.getDistance());
			String nname=n.getId().trim();
			double[] cn=(double[])hm.get(nname);
			String c1name=n.getLeft().getId().trim();
			double[] c1=(double[])hm.get(c1name);
			String c2name=n.getRight().getId().trim();
			double[] c2=(double[])hm.get(c2name);
			double[] X=new double[3];
			X[0]=-1;
			X[1]=this.dis(cn, c1);
			X[2]=this.dis(cn, c2);
			hm_mds_innernodes_dis.put(n.getId().trim(), X);
			//System.out.println("MDSDis:"+X[0]+","+X[1]+","+X[2]);
			node lnode=n.getLeft();
			if(lnode.getLeft()!=null && lnode.getRight()!=null){
				String lname=lnode.getId().trim();
				double[] cL=(double[])hm.get(lname);
				double[] XL=(double[])hm_mds_innernodes_dis.get(lname);
				XL[0]=this.dis(cL, cn);
				//System.out.println("MDSDis:"+XL[0]);
				hm_mds_innernodes_dis.put(lname, XL);
			}
			
			node rnode=n.getRight();
			if(rnode.getLeft()!=null && rnode.getRight()!=null){
				String rname=rnode.getId().trim();
				double[] cR=(double[])hm.get(rname);
				double[] XR=(double[])hm_mds_innernodes_dis.get(rname);
				XR[0]=this.dis(cR, cn);
				//System.out.println("MDSDis:"+XR[0]);
				hm_mds_innernodes_dis.put(rname, XR);
			}
			//System.out.println("MDSDis:"+X[0]+","+X[1]+","+X[2]);
			//X=(double[])hm_mds_innernodes_dis.get(nname);
			//System.out.println("MDSDis:"+X[0]+","+X[1]+","+X[2]);
		}else{
			//System.out.println(n.getId()+":"+n.getDistance());
			return;
		}
		
		
	}
	
	
	//In use 
	//Generate distance matrix for shannon mapping 
	public void postTree_cal_mds_dis2(node n, HashMap hm){
		if(n.getLeft()!=null && n.getRight()!=null){
			postTree_cal_mds_dis2(n.getLeft(),hm);
			postTree_cal_mds_dis2(n.getRight(),hm);
			String nname=n.getId().trim();
			double[] cn=(double[])hm.get(nname);
			String c1name=n.getLeft().getId().trim();
			double[] c1=(double[])hm.get(c1name);
			String c2name=n.getRight().getId().trim();
			double[] c2=(double[])hm.get(c2name);
			double[] X=new double[3];
			X[1]=this.dis(cn, c1);
			X[2]=this.dis(cn, c2);
			String namec=n.getId().trim();
			String namel=n.getLeft().getId().trim();
			String namer=n.getRight().getId().trim();
			hm_mds_innernodes_dis2.put(namec+":"+namel, new Double(X[1]));
			hm_mds_innernodes_dis2.put(namel+":"+namec, new Double(X[1]));
			hm_mds_innernodes_dis2.put(namec+":"+namer, new Double(X[2]));
			hm_mds_innernodes_dis2.put(namer+":"+namec, new Double(X[2]));
		}else{
			return;
		}
	}
	
	
	//In use 
	//Generate distance maxtrix - Hashmap for shannon mapping, will call postTree_cal_mds_dis2 and add root
	public void cal_mds_dis(node root, HashMap name_coords){
		this.hm_mds_innernodes_dis2=new HashMap();
		this.postTree_cal_mds_dis2(root.getLeft(), name_coords);
		String cname=root.getId().trim();
		String lname=root.getLeft().getId().trim();
		double[] cRoot=(double[])name_coords.get(cname);
		double[] cLeft=(double[])name_coords.get(lname);
		//todo: potential bugs
		if(cRoot==null || cLeft==null){
			double XL= 0;
			this.hm_mds_innernodes_dis2.put(cname+":"+lname, new Double(XL));
			this.hm_mds_innernodes_dis2.put(lname+":"+cname, new Double(XL));
		}else{
			double XL=this.dis(cRoot, cLeft);
			this.hm_mds_innernodes_dis2.put(cname+":"+lname, new Double(XL));
			this.hm_mds_innernodes_dis2.put(lname+":"+cname, new Double(XL));

		}
	}
	
	public void postTreefindConnetedNodeName(node n){
		if(n.getLeft()!=null && n.getRight()!=null){
			postTreefindConnetedNodeName(n.getLeft());
			postTreefindConnetedNodeName(n.getRight());
			//System.out.println(n.getId()+":"+n.getDistance());
			String[] names=new String[3];
			names[1]=n.getLeft().getId().trim();
			names[2]=n.getRight().getId().trim();
			
			hm_node_connectednodename.put(n.getId().trim(), names);
			
			node lnode=n.getLeft();
			if(lnode.getLeft()!=null && lnode.getRight()!=null){
				String lname=lnode.getId().trim();
				
				String[] namesL=(String[])hm_node_connectednodename.get(lname);
				
				namesL[0]=n.getId().trim();
				
				hm_node_connectednodename.put(lname, namesL);
			}
			
			node rnode=n.getRight();
			if(rnode.getLeft()!=null && rnode.getRight()!=null){
				String rname=rnode.getId().trim();
				
				String[] namesR=(String[])hm_node_connectednodename.get(rname);
				
				namesR[0]=n.getId().trim();
				
				hm_node_connectednodename.put(rname, namesR);
			}
			
			
		}else{
			//System.out.println(n.getId()+":"+n.getDistance());
			return;
		}
		
	}
	
	public HashMap getTreeNameToDis(){
		return this.hm_innernodes_dis;
	}
	
	public HashMap getMDSNameToDis(){
		return this.hm_mds_innernodes_dis;
	}
	
	public HashMap getTreeNameToDis2(){
		return this.hm_innernodes_dis2;
	}
	
	public HashMap getMDSNameToDis2(){
		return this.hm_mds_innernodes_dis2;
	}
	
	public HashMap getNodeToConnetednodename(){
		return this.hm_node_connectednodename;
	}
	
	public HashMap getControlp(){
		return this.hm_controlpoint;
	}
	
	public void printHm(){
		Set s=hm_innernodes_dis.keySet();
		Iterator it=s.iterator();
		while(it.hasNext()){
			String key=(String)it.next();
			double[] X=(double[])hm_innernodes_dis.get(key);
			System.out.println(key+":"+X[0]+" "+X[1]+" "+X[2]);
		}
	}
	
	public void printHm(HashMap hm){
		Set s=hm.keySet();
		Iterator it=s.iterator();
		while(it.hasNext()){
			String key=(String)it.next();
			double[] X=(double[])hm.get(key);
			System.out.println(key+":"+X[0]+" "+X[1]+" "+X[2]);
		}
	}
	
	public void printHm2(HashMap hm){
		Set s=hm.keySet();
		Iterator it=s.iterator();
		while(it.hasNext()){
			String key=(String)it.next();
			String[] X=(String[])hm.get(key);
			System.out.println(key+":"+X[0]+" "+X[1]+" "+X[2]);
		}
	}
	
	public void printHm3(HashMap hm){
		Set s=hm.keySet();
		Iterator it=s.iterator();
		while(it.hasNext()){
			String key=(String)it.next();
			Double X=(Double)hm.get(key);
			System.out.println(key+":"+X);
		}
	}
	
	public void genDisMatrixFromTree(node n, File fout) throws Exception{
		this.postTree(n.getLeft());
		this.nodes_list.add(n.getId());
		this.postTree_cal_tree_dis2(n.getLeft());
		String namec=n.getId().trim();
		String namel=n.getLeft().getId().trim();
		this.hm_innernodes_dis2.put(namec+":"+namel, new Double(n.getLeft().getDistance()));
		this.hm_innernodes_dis2.put(namel+":"+namec, new Double(n.getLeft().getDistance()));
		
		StringBuffer sb=new StringBuffer("");
		for(int i=0;i<this.nodes_list.size();i++){
			String node_name=(String)this.nodes_list.get(i);
			for(int k=0;k<this.nodes_list.size();k++){
				String temp=(String)this.nodes_list.get(k);
				if(i!=k){
					if(hm_innernodes_dis2.containsKey(node_name+":"+temp)){
						Double dd=(Double)hm_innernodes_dis2.get(node_name+":"+temp);
						sb.append(dd.doubleValue()+"	");
						System.out.print("-,");
					}else if(hm_innernodes_dis2.containsKey(temp+":"+node_name)){
						Double dd=(Double)hm_innernodes_dis2.get(node_name+":"+temp);
						sb.append(dd.doubleValue()+"	");
						System.out.print("-,");
					}else{
						sb.append("NaN"+"	");
					}
				}else{
					
					sb.append("0"+"	");
				}
			}
			System.out.println();
			sb.append("\n");
		}
		
		FileWriter fw=new FileWriter(fout);
		fw.write(sb.toString());
		fw.close();
		
	}
	
	
	//In use 
	//Fill the node_list and tree distances
	//return sum of all tree distances
	public double sumTreeSpacedistance(node n){
		
		//fill in node_list and tree distance matrix
		double c=0;
		this.postTree(n.getLeft());
		this.nodes_list.add(n.getId());
		
		this.postTree_cal_tree_dis2(n.getLeft());
		String namec=n.getId().trim();
		String namel=n.getLeft().getId().trim();
		this.hm_innernodes_dis2.put(namec+":"+namel, new Double(n.getLeft().getDistance()*scaleFactor));
		this.hm_innernodes_dis2.put(namel+":"+namec, new Double(n.getLeft().getDistance()*scaleFactor));
		
		for(int i=0;i<this.nodes_list.size();i++){
			String node_name=(String)this.nodes_list.get(i);
			for(int k=0;k<this.nodes_list.size();k++){
				String temp=(String)this.nodes_list.get(k);
					if(hm_innernodes_dis2.containsKey(node_name+":"+temp)){
						Double dd=(Double)hm_innernodes_dis2.get(node_name+":"+temp);
						c=c+dd.doubleValue();
					}
			}
			
		}
		return c;
	}
	
	
	//In use 
	//Calculate error
	public double calError(double c){
		double err=0;
		for(int i=0;i<this.nodes_list.size();i++){
			String node_name=(String)this.nodes_list.get(i);
			for(int k=0;k<this.nodes_list.size();k++){
				String temp=(String)this.nodes_list.get(k);
					if(hm_innernodes_dis2.containsKey(node_name+":"+temp)){
						Double D_tree=(Double)hm_innernodes_dis2.get(node_name+":"+temp);
						Double D_mds=(Double)hm_mds_innernodes_dis2.get(node_name+":"+temp);
						double d_tree=D_tree.doubleValue();
						double d_mds=D_mds.doubleValue();
						err=err+(d_tree-d_mds)*(d_tree-d_mds)/d_tree;
						if(d_tree<=0){
							System.out.println(node_name+":"+temp+":"+d_tree);
						}
					}
			}
			
		}
		return err/c;
	}
	
	public double calErrorBezier(double c){
		//this.printHm3(hm_innernodes_dis2);
		//System.out.println("-----------");
		//this.printHm3(hm_mds_innernodes_dis2);
		
		double err=0;
		for(int i=0;i<this.nodes_list.size();i++){
			String node_name=(String)this.nodes_list.get(i);
			for(int k=i;k<this.nodes_list.size();k++){
				String temp=(String)this.nodes_list.get(k);
				//if(i<k){
					if(hm_innernodes_dis2.containsKey(node_name+":"+temp)){
						Double D_tree=(Double)hm_innernodes_dis2.get(node_name+":"+temp);
						Double D_mds=(Double)hm_mds_innernodes_dis2.get(node_name+":"+temp);
						
						double[] D_bezier=(double[])this.hm_controlpoint.get(node_name+":"+temp);
						if(D_bezier!=null){
							double d_bezier=D_bezier[2];
						
							double d_tree=D_tree.doubleValue();
							double d_mds=d_bezier;
							err=err+(d_tree-d_mds)*(d_tree-d_mds)/d_tree;
						}
					}
				//}
			}
			
		}
		
		return err/c;
	}
	
	public void ccPlot(File fout) throws Exception{
		
		StringBuffer sb=new StringBuffer("");
		for(int i=0;i<this.nodes_list.size();i++){
			String node_name=(String)this.nodes_list.get(i);
			for(int k=i;k<this.nodes_list.size();k++){
				String temp=(String)this.nodes_list.get(k);
				//if(i<k){
					if(hm_innernodes_dis2.containsKey(node_name+":"+temp)){
						Double D_tree=(Double)hm_innernodes_dis2.get(node_name+":"+temp);
						Double D_mds=(Double)hm_mds_innernodes_dis2.get(node_name+":"+temp);
						double d_tree=D_tree.doubleValue();
						double d_mds=D_mds.doubleValue();
						sb.append(d_tree+"	"+d_mds+"\n");
						
					}
				//}
			}
			
		}
		
		FileWriter fw=new FileWriter(fout);
		fw.write(sb.toString());
		fw.close();
	}
	
	public void ccPlotBezier(File fout) throws Exception{
		
		StringBuffer sb=new StringBuffer("");
		for(int i=0;i<this.nodes_list.size();i++){
			String node_name=(String)this.nodes_list.get(i);
			for(int k=i;k<this.nodes_list.size();k++){
				String temp=(String)this.nodes_list.get(k);
				//if(i<k){
					if(hm_innernodes_dis2.containsKey(node_name+":"+temp)){
						Double D_tree=(Double)hm_innernodes_dis2.get(node_name+":"+temp);
						double[] D_bezier=(double[])this.hm_controlpoint.get(node_name+":"+temp);
						if(D_tree!=null && D_bezier!=null){
							double d_tree=D_tree.doubleValue();
							double d_bezier=D_bezier[2];
							sb.append(d_tree+"	"+d_bezier+"\n");
						}
					}
				//}
			}
			
		}
		
		FileWriter fw=new FileWriter(fout);
		fw.write(sb.toString());
		fw.close();
	}
	
	public ArrayList ccMatrix() throws Exception{
		
		ArrayList al=new ArrayList();
		
		for(int i=0;i<this.nodes_list.size();i++){
			String node_name=(String)this.nodes_list.get(i);
			for(int k=i;k<this.nodes_list.size();k++){
				String temp=(String)this.nodes_list.get(k);
				//if(i<k){
					if(hm_innernodes_dis2.containsKey(node_name+":"+temp)){
						Double D_tree=(Double)hm_innernodes_dis2.get(node_name+":"+temp);
						Double D_mds=(Double)hm_mds_innernodes_dis2.get(node_name+":"+temp);
						double d_tree=D_tree.doubleValue();
						double d_mds=D_mds.doubleValue();
						double[] D=new double[2];
						D[0]=d_tree;
						D[1]=d_mds;
						al.add(D);	
						
					}
				//} v 
			}
			
		}
		
		return al;
	}
	
	public ArrayList ccMatrixBezier() throws Exception{
		
		ArrayList al=new ArrayList();
		
		for(int i=0;i<this.nodes_list.size();i++){
			String node_name=(String)this.nodes_list.get(i);
			for(int k=i;k<this.nodes_list.size();k++){
				String temp=(String)this.nodes_list.get(k);
				//if(i<k){
					if(hm_innernodes_dis2.containsKey(node_name+":"+temp)){
						Double D_tree=(Double)hm_innernodes_dis2.get(node_name+":"+temp);
						double d_tree=D_tree.doubleValue();
						
						double[] D_bezier=(double[])this.hm_controlpoint.get(node_name+":"+temp);
						if(D_bezier!=null){
							double d_bezier=D_bezier[2];
							
							double[] D=new double[2];
							D[0]=d_tree;
							D[1]=d_bezier;
							al.add(D);
						}
					}
				//} v 
			}
			
		}
		
		return al;
	}
	
	public double Calcc(ArrayList al){
		double cc=0;
		int n=al.size();
		double[] x=new double[n];
		double[] y=new double[n];
		double[] xy=new double[n];
		double[] x2=new double[n];
		double[] y2=new double[n];
		for(int i=0;i<n;i++){
			double[] temp=(double[])al.get(i);
			x[i]=temp[0];
			y[i]=temp[1];
			xy[i]=temp[0]*temp[1];
			x2[i]=temp[0]*temp[0];
			y2[i]=temp[1]*temp[1];
		}
		double sx=this.sum(x);
		double sy=this.sum(y);
		double sxy=this.sum(xy);
		double sx2=this.sum(x2);
		double sy2=this.sum(y2);
		
		cc=(n*sxy-sx*sy)/Math.sqrt((n*sx2-sx*sx)*(n*sy2-sy*sy));
		
		return cc;
	}
	
	public double rmsPositive(ArrayList al, double slop){
		double rms=0;
		int cnt=0;
		for(int i=0;i<al.size();i++){
			double[] temp=(double[])al.get(i);
			if(temp[1]>temp[0]*slop){
				rms=rms+(temp[1]-temp[0]*slop)*(temp[1]-temp[0]*slop);
				cnt++;
			}
		}
		
		
		return rms/(double)cnt;
	}
	
	public double[] linearRegression(ArrayList al){
		
		double[] parameter=new double[2]; // y=mx+b
		double m=0;
		double b=0;
		int n=al.size();
		double[] x=new double[n];
		double[] y=new double[n];
		double[] xy=new double[n];
		double[] x2=new double[n];
		double[] y2=new double[n];
		for(int i=0;i<n;i++){
			double[] temp=(double[])al.get(i);
			x[i]=temp[0];
			y[i]=temp[1];
			xy[i]=temp[0]*temp[1];
			x2[i]=temp[0]*temp[0];
			y2[i]=temp[1]*temp[1];
		}
		double sx=this.sum(x);
		double sy=this.sum(y);
		double sxy=this.sum(xy);
		double sx2=this.sum(x2);
		double sy2=this.sum(y2);                              
		
		m=(n*sxy-sy*sx)/(n*sx2-sx*sx);
		b=(sy-m*sx)/n;
		
		parameter[0]=m;
		parameter[1]=b;
		
		return parameter;
	}
	
	public double[] multiLinearRegression(ArrayList al, int n){
		double[] parameter=new double[2]; // y=mx+b
		ArrayList alt=al;
		parameter=this.linearRegression(alt);
		for(int i=0;i<(n-1);i++){
			
			for(int k=0;k<alt.size();k++){
				double[] temp=(double[])alt.get(k);
				double ty=temp[0]*parameter[0];
				if(temp[1]<ty){
					temp[0]=temp[0]*parameter[0];
				}
				alt.set(k, temp);
			}
			parameter=this.linearRegression(alt);
			//System.out.println(parameter[0]);
		}
		return parameter;
	}
	
	public double adaline(ArrayList al){
		double w=0;
		int k_max=100;
		double ebsilon_s=0.1;
		double ebsilon_e=0.001;
		double ebsilon=ebsilon_s;
		double y=0;
		double dw=0;
		
		for(int i=0;i<k_max;i++){
			
			for(int k=0;k<al.size();k++){
				double[] temp=(double[])al.get(k);
				double x=temp[0];
				double s=temp[1];
				y=w*x;
				dw=ebsilon*(s-y)*x;
				w+=dw;
				ebsilon=ebsilon_s*Math.pow((ebsilon_e/ebsilon_s), (i/k_max));
			}
			
		}
		
		return w;
	}
	
	public double multiAdaline(ArrayList al, int n){
		
		ArrayList alt=al;
		double parameter=this.adaline(alt);
		for(int i=0;i<(n-1);i++){
			
			for(int k=0;k<alt.size();k++){
				double[] temp=(double[])alt.get(k);
				double ty=temp[0]*parameter;
				if(temp[1]<ty){
					temp[0]=temp[0]*parameter;
				}
				alt.set(k, temp);
			}
			parameter=this.adaline(alt);
			//System.out.println(parameter);
		}
		return parameter;
	}
	
	public double sum(double[] dl){
		double d=0;
		for(int i=0;i<dl.length;i++){
			d+=dl[i];
		}
		return d;
	}
	
	public HashMap upDate(String nodeName, HashMap Name_Coords){
		
		pnp Pnp=new pnp();
		double d1x=0;
		double d1y=0;
		double d2x=0;
		double d2y=0;
		double[] coord_upnode=(double[])Name_Coords.get(nodeName);
		for(int i=0;i<this.nodes_list.size();i++){
			String cnode=(String)this.nodes_list.get(i);
			Double D=(Double)this.hm_innernodes_dis2.get(nodeName+":"+cnode);
			if(D!=null){
			
				double d_tree=((Double)this.hm_innernodes_dis2.get(nodeName+":"+cnode)).doubleValue();
				double d_mds=((Double)this.hm_mds_innernodes_dis2.get(nodeName+":"+cnode)).doubleValue();
				double[] coord_cnode=(double[])Name_Coords.get(cnode);
				if (coord_cnode != null){
					d1x=d1x+this.firstDstep(d_tree, d_mds, coord_cnode[0], coord_upnode[0]);
					//d2x=d2x+this.firstDstep(d_tree, d_mds, coord_cnode[0], coord_upnode[0]);
					
					d2x=d2x+this.secondDstep(d_tree, d_mds, coord_cnode[0], coord_upnode[0]);
					
					d1y=d1y+this.firstDstep(d_tree, d_mds, coord_cnode[1], coord_upnode[1]);
					//d2y=d2y+this.firstDstep(d_tree, d_mds, coord_cnode[1], coord_upnode[1]);
					
					d2y=d2y+this.secondDstep(d_tree, d_mds, coord_cnode[1], coord_upnode[1]);
				}
			}
		}
		if(d2x<0){
			d2x=-d2x;	
		}
		
		if(d2y<0){
			d2y=-d2y;
		}
		
		double deltaX=-d1x/d2x;
		
		coord_upnode[0]=coord_upnode[0]-mf*deltaX;
		double deltaY=-d1y/d2y;
		
		coord_upnode[1]=coord_upnode[1]-mf*deltaY;	
		
		Name_Coords.put(nodeName, coord_upnode);
		return Name_Coords;
	}
	
	
	public HashMap upDateHD(String nodeName, HashMap Name_Coords, int dim){
		
		pnp Pnp=new pnp();
		double[] d1x=new double[dim];
		double[] d2x=new double[dim];
		for(int i=0;i<dim;i++){
			d1x[i]=0;
			d2x[i]=0;
		}
		//double d1x=0;
		double d1y=0;
		//double d2x=0;
		double d2y=0;
		double[] coord_upnode=(double[])Name_Coords.get(nodeName);
		for(int i=0;i<this.nodes_list.size();i++){
			String cnode=(String)this.nodes_list.get(i);
			Double D=(Double)this.hm_innernodes_dis2.get(nodeName+":"+cnode);
			if(D!=null){
			
				double d_tree=((Double)this.hm_innernodes_dis2.get(nodeName+":"+cnode)).doubleValue();
				double d_mds=((Double)this.hm_mds_innernodes_dis2.get(nodeName+":"+cnode)).doubleValue();
				double[] coord_cnode=(double[])Name_Coords.get(cnode);
				
				for(int d=0;d<coord_cnode.length;d++){
					d1x[d]=d1x[d]+this.firstDstep(d_tree, d_mds, coord_cnode[d], coord_upnode[d]);
					//d2x=d2x+this.firstDstep(d_tree, d_mds, coord_cnode[0], coord_upnode[0]);
					
					d2x[d]=d2x[d]+this.secondDstep(d_tree, d_mds, coord_cnode[d], coord_upnode[d]);
				}
				
				//d1x=d1x+this.firstDstep(d_tree, d_mds, coord_cnode[0], coord_upnode[0]);
				//d2x=d2x+this.firstDstep(d_tree, d_mds, coord_cnode[0], coord_upnode[0]);
				
				//d2x=d2x+this.secondDstep(d_tree, d_mds, coord_cnode[0], coord_upnode[0]);
				
				//d1y=d1y+this.firstDstep(d_tree, d_mds, coord_cnode[1], coord_upnode[1]);
				//d2y=d2y+this.firstDstep(d_tree, d_mds, coord_cnode[1], coord_upnode[1]);
				
				//d2y=d2y+this.secondDstep(d_tree, d_mds, coord_cnode[1], coord_upnode[1]);
			
			}
		}
		
		for(int d=0;d<dim;d++){
			if(d2x[d]<0){
				d2x[d]=-d2x[d];	
			}
			double deltaX=-d1x[d]/d2x[d];
			coord_upnode[d]=coord_upnode[d]-mf*deltaX;
			
		
		}
		
		//if(d2y<0){
		//	d2y=-d2y;
		//}
		
		//double deltaX=-d1x/d2x;
		
		//coord_upnode[0]=coord_upnode[0]-mf*deltaX;
		//double deltaY=-d1y/d2y;
		
		//coord_upnode[1]=coord_upnode[1]-mf*deltaY;	
		
		Name_Coords.put(nodeName, coord_upnode);
		return Name_Coords;
	}
	
	public HashMap upDateBelow(String nodeName, HashMap Name_Coords){
		
		pnp Pnp=new pnp();
		double d1x=0;
		double d1y=0;
		double d2x=0;
		double d2y=0;
		double[] coord_upnode=(double[])Name_Coords.get(nodeName);
		boolean flag=false;
		for(int i=0;i<this.nodes_list.size();i++){
			String cnode=(String)this.nodes_list.get(i);
			Double D=(Double)this.hm_innernodes_dis2.get(nodeName+":"+cnode);
			if(D!=null){
			
				double d_tree=((Double)this.hm_innernodes_dis2.get(nodeName+":"+cnode)).doubleValue();
				double d_mds=((Double)this.hm_mds_innernodes_dis2.get(nodeName+":"+cnode)).doubleValue();
				double[] coord_cnode=(double[])Name_Coords.get(cnode);
				
				if(d_tree<d_mds){
				flag=true;
				d1x=d1x+this.firstDstep(d_tree, d_mds, coord_cnode[0], coord_upnode[0]);
				//d2x=d2x+this.firstDstep(d_tree, d_mds, coord_cnode[0], coord_upnode[0]);
				d2x=d2x+this.secondDstep(d_tree, d_mds, coord_cnode[0], coord_upnode[0]);
				
				
				d1y=d1y+this.firstDstep(d_tree, d_mds, coord_cnode[1], coord_upnode[1]);
				//d2y=d2y+this.firstDstep(d_tree, d_mds, coord_cnode[1], coord_upnode[1]);
				d2y=d2y+this.secondDstep(d_tree, d_mds, coord_cnode[1], coord_upnode[1]);
				
				}
			}
		}
		if(d2x<0){
			d2x=-d2x;	
		}
		
		if(d2y<0){
			d2y=-d2y;
		}
		
		if(flag){
			double deltaX=-d1x/d2x;
			
			coord_upnode[0]=coord_upnode[0]-mf*deltaX;
			double deltaY=-d1y/d2y;
			
			coord_upnode[1]=coord_upnode[1]-mf*deltaY;	
			
			Name_Coords.put(nodeName, coord_upnode);
		}
		return Name_Coords;
	}
	
	
	public HashMap upDateBelowHD(String nodeName, HashMap Name_Coords, int dim){
		
		pnp Pnp=new pnp();
		double[] d1x=new double[dim];
		double[] d2x=new double[dim];
		for(int i=0;i<dim;i++){
			d1x[i]=0;
			d2x[i]=0;
		}
		
		double d1y=0;
		double d2y=0;
		double[] coord_upnode=(double[])Name_Coords.get(nodeName);
		boolean flag=false;
		for(int i=0;i<this.nodes_list.size();i++){
			String cnode=(String)this.nodes_list.get(i);
			Double D=(Double)this.hm_innernodes_dis2.get(nodeName+":"+cnode);
			if(D!=null){
			
				double d_tree=((Double)this.hm_innernodes_dis2.get(nodeName+":"+cnode)).doubleValue();
				double d_mds=((Double)this.hm_mds_innernodes_dis2.get(nodeName+":"+cnode)).doubleValue();
				double[] coord_cnode=(double[])Name_Coords.get(cnode);
				
				if(d_tree<d_mds){
					flag=true;
				
					for(int d=0;d<coord_cnode.length;d++){
						d1x[d]=d1x[d]+this.firstDstep(d_tree, d_mds, coord_cnode[d], coord_upnode[d]);
						
						
						d2x[d]=d2x[d]+this.secondDstep(d_tree, d_mds, coord_cnode[d], coord_upnode[d]);
					}
				}

			
			}
		}
		
		if(flag){
			for(int d=0;d<dim;d++){
				if(d2x[d]<0){
					d2x[d]=-d2x[d];	
				}
				double deltaX=-d1x[d]/d2x[d];
				coord_upnode[d]=coord_upnode[d]-mf*deltaX;
			}
			Name_Coords.put(nodeName, coord_upnode);
		}
		return Name_Coords;
	}
	
	
	public HashMap upDateSmall(String nodeName, HashMap Name_Coords){
		
		pnp Pnp=new pnp();
		double d1x=0;
		double d1y=0;
		double d2x=0;
		double d2y=0;
		double[] coord_upnode=(double[])Name_Coords.get(nodeName);
		for(int i=0;i<this.nodes_list.size();i++){
			String cnode=(String)this.nodes_list.get(i);
			Double D=(Double)this.hm_innernodes_dis2.get(nodeName+":"+cnode);
			if(D!=null){
			
				double d_tree=((Double)this.hm_innernodes_dis2.get(nodeName+":"+cnode)).doubleValue();
				double d_mds=((Double)this.hm_mds_innernodes_dis2.get(nodeName+":"+cnode)).doubleValue();
				double[] coord_cnode=(double[])Name_Coords.get(cnode);
				
				d1x=d1x+this.firstDstepSmall(d_tree, d_mds, coord_cnode[0], coord_upnode[0]);
				d2x=d2x+this.firstDstepSmall(d_tree, d_mds, coord_cnode[0], coord_upnode[0]);
				
				d1y=d1y+this.firstDstepSmall(d_tree, d_mds, coord_cnode[1], coord_upnode[1]);
				d2y=d2y+this.firstDstepSmall(d_tree, d_mds, coord_cnode[1], coord_upnode[1]);
				
			
			}
		}
		if(d2x<0){
			d2x=-d2x;	
		}
		
		if(d2y<0){
			d2y=-d2y;
		}
		
		double deltaX=-d1x/d2x;
		
		coord_upnode[0]=coord_upnode[0]-mf*deltaX;
		double deltaY=-d1y/d2y;
		
		coord_upnode[1]=coord_upnode[1]-mf*deltaY;	
		
		Name_Coords.put(nodeName, coord_upnode);
		return Name_Coords;
	}
	
	public double firstDstep(double d_tree, double d_mds, double d_coord_c, double d_coord_u){
		double d=0;
		d=((d_tree-d_mds)/(d_tree*d_mds))*(d_coord_u-d_coord_c);
		return d;
	}
	
	public double firstDstepSmall(double d_tree, double d_mds, double d_coord_c, double d_coord_u){
		double d=0;
		if(d_tree-d_mds<0){
			double penalty=(d_mds-d_tree)*10000;
			d=penalty*((d_tree-d_mds)/(d_tree*d_mds))*(d_coord_u-d_coord_c);
		}else{
			d=((d_tree-d_mds)/(d_tree*d_mds))*(d_coord_u-d_coord_c);
		}
		return d;
	}
	
	public double secondDstepSmall(double d_tree, double d_mds, double d_coord_c, double d_coord_u){
		double d=0;
		if(d_tree-d_mds<0){
			double penalty=(d_mds-d_tree)*10000;
			d=penalty*((d_tree-d_mds)-((d_coord_u-d_coord_c)*(d_coord_u-d_coord_c)/d_mds)*(1+((d_tree-d_mds)/d_mds)))/(d_tree*d_mds);
		}else{
			d=((d_tree-d_mds)-((d_coord_u-d_coord_c)*(d_coord_u-d_coord_c)/d_mds)*(1+((d_tree-d_mds)/d_mds)))/(d_tree*d_mds);
		}
		//d=((d_tree-d_mds)-((d_coord_u-d_coord_c)*(d_coord_u-d_coord_c)/d_mds)*(1+((d_tree-d_mds)/d_mds)))/(d_tree*d_mds);
		return d;
	}
	
	public double secondDstep(double d_tree, double d_mds, double d_coord_c, double d_coord_u){
		double d=0;
		d=((d_tree-d_mds)-((d_coord_u-d_coord_c)*(d_coord_u-d_coord_c)/d_mds)*(1+((d_tree-d_mds)/d_mds)))/(d_tree*d_mds);
		return d;
	}
	
	public void genTreePlot(File fout, node root, HashMap Name_Coords,HashMap hm_controlpoint) throws Exception{
		
		
		this.postPlot(root.getLeft(), Name_Coords,hm_controlpoint);
		String namec=root.getId().trim();
		String namel=root.getLeft().getId().trim();
		double[] Dc=(double[])Name_Coords.get(namec);
		double[] Dl=(double[])Name_Coords.get(namel);
		double[] CPl=(double[])hm_controlpoint.get(namec+":"+namel);
		if(Dc!=null && Dl!=null && CPl!=null){
			sb.append(Dc[0]+"	"+Dc[1]+"	"+Dl[0]+"	"+Dl[1]+"	"+CPl[0]+"	"+CPl[1]+"\n");
		}
		FileWriter fw=new FileWriter(fout);
		fw.write(sb.toString());
		fw.close();
		
		
	}
	
	public void postPlot(node n, HashMap Name_Coords, HashMap hm_controlpoint){
		if(n.getLeft()!=null && n.getRight()!=null){
			postPlot(n.getLeft(),Name_Coords,hm_controlpoint);
			postPlot(n.getRight(),Name_Coords,hm_controlpoint);
			String namec=n.getId().trim();
			String namel=n.getLeft().getId().trim();
			String namer=n.getRight().getId().trim();
			
			double[] Dc=(double[])Name_Coords.get(namec);
			double[] Dl=(double[])Name_Coords.get(namel);
			double[] Dr=(double[])Name_Coords.get(namer);
			double[] CPl=(double[])hm_controlpoint.get(namec+":"+namel);
			double[] CPr=(double[])hm_controlpoint.get(namec+":"+namer);
			sb.append(Dc[0]+"	"+Dc[1]+"	"+Dl[0]+"	"+Dl[1]+"	"+CPl[0]+"	"+CPl[1]+"\n");
			sb.append(Dc[0]+"	"+Dc[1]+"	"+Dr[0]+"	"+Dr[1]+"	"+CPr[0]+"	"+CPr[1]+"\n");
			
			
			
		}else{
			
			return;
		}
		
	}
	
	public void postTreeBezier(node n, double scale, HashMap Name_Coords){
		if(n.getLeft()!=null && n.getRight()!=null){
			this.postTreeBezier(n.getLeft(), scale, Name_Coords);
			this.postTreeBezier(n.getRight(), scale, Name_Coords);
			String nname=n.getId().trim();
			double[] cn=(double[])Name_Coords.get(nname);
			String c1name=n.getLeft().getId().trim();
			double[] c1=(double[])Name_Coords.get(c1name);
			String c2name=n.getRight().getId().trim();
			double[] c2=(double[])Name_Coords.get(c2name);
			
			Double Dtl=(Double)this.hm_innernodes_dis2.get(nname+":"+c1name);
			Double Dtr=(Double)this.hm_innernodes_dis2.get(nname+":"+c2name);
			
			//double[] contrlPL=PNP.findControlPoint(cn[0], cn[1], c1[0], c1[1], scale, Dtl.doubleValue(), error);
			//double[] contrlPR=PNP.findControlPoint(cn[0], cn[1], c2[0], c2[1], scale, Dtr.doubleValue(), error);
			
			double[] contrlPL=vis.segSpace(cn[0], cn[1], c1[0], c1[1], scale, Dtl.doubleValue(), error);
			double[] contrlPR=vis.segSpace(cn[0], cn[1], c2[0], c2[1], scale, Dtr.doubleValue(), error);
			
			hm_controlpoint.put(nname+":"+c1name, contrlPL);
			hm_controlpoint.put(nname+":"+c2name, contrlPR);
			
			hm_controlpoint.put(c1name+":"+nname, contrlPL);
			hm_controlpoint.put(c2name+":"+nname, contrlPR);
		}else{
			return;
		}
	}
	
	public void calBezier(node root, HashMap Name_Coords, double scale){
		//double scale=PNP.maxScale(this.hm_innernodes_dis2, this.hm_mds_innernodes_dis2);
		//System.out.println(scale);
		HashMap N_Name_Coords=(HashMap)Name_Coords.clone();
		vis.init(N_Name_Coords);
		this.postTreeBezier(root.getLeft(), scale, Name_Coords);
		String cname=root.getId().trim();
		String lname=root.getLeft().getId().trim();
		double[] cn=(double[])Name_Coords.get(cname);
		double[] c1=(double[])Name_Coords.get(lname);
		
		if(cn!=null && c1!=null){
			Double Dtl=(Double)this.hm_innernodes_dis2.get(cname+":"+lname);
			//double[] contrlPL=PNP.findControlPoint(cn[0], cn[1], c1[0], c1[1], scale, Dtl.doubleValue(), error);
			double[] contrlPL=vis.segSpace(cn[0], cn[1], c1[0], c1[1], scale, Dtl.doubleValue(), error);
			hm_controlpoint.put(cname+":"+lname, contrlPL);
			hm_controlpoint.put(lname+":"+cname, contrlPL);
		}
	}
	
	public static void main(String[] args) {
		// TODO Auto-generated method stub

		
	}

}

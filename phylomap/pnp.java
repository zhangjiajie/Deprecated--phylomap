package phylomap;
import java.util.*;
import java.io.*;
import db.*;


public class pnp {

	/**
	 * @param args
	 */
	double mf=0.01;
	public node tree_root=null;
	public HashMap hm_controlpoint=new HashMap();
	public double sammon_error=0;
	public double bezier_error=0;
	public double cc=0;
	public double bezier_cc=0;
	
	//In use
	//Random assign coordinates to inner nodes of the tree to double[][]
	public double[][] init(int numPoint, double xlim1, double xlim2, double ylim1, double ylim2){
		Date date=new Date();
		Random rand=new Random(date.getTime());
		double[][] coords=new double[numPoint][2];
		for(int i=0;i<numPoint;i++){
			double dx=xlim1+(xlim2-xlim1)*rand.nextDouble();
			double dy=ylim1+(ylim2-ylim1)*rand.nextDouble();
			
			coords[i][0]=dx;
			coords[i][1]=dy;
		}
		return coords;
	}
	
	//In use 
	//Load leaf coordinates from the PCoA result to HashMap
	public HashMap loadLeafCoords(File fin) throws Exception{
		HashMap hm=new HashMap();
		FileReader fr=new FileReader(fin);
		BufferedReader bf=new BufferedReader(fr);
		String s=bf.readLine();
		while(s!=null){
			s=s.trim();
			
			String[] ss=s.split("	");
			
			double dx=Double.parseDouble(ss[1]);
			double dy=Double.parseDouble(ss[2]);
			double[] da=new double[2];
			da[0]=dx;
			da[1]=dy;
			hm.put(ss[0].trim(), da);
			s=bf.readLine();
		}
		
		bf.close();
		fr.close();
		return hm;
		
	}
	
	//In use 
	//Assign inti coords to inner nodes, return hashmap
	public HashMap assignCoordsToInnernodes(HashMap Leaf_Name_Coords, double[][] coords){
		HashMap nodes_coords=new HashMap();
		Set s=Leaf_Name_Coords.keySet();
		Iterator it=s.iterator();
		int cnt=0;
		while(it.hasNext()){
			nodes_coords.put(it.next(), coords[cnt++]);
		}
		return nodes_coords;
	}
	
	//In use 
	//Generate a bigger HashMap from 2 samll hashmap
	public HashMap intergateHMs(HashMap h1, HashMap h2){
		HashMap h3=new HashMap();
		Set s=h1.keySet();
		Iterator it=s.iterator();
		while(it.hasNext()){
			Object o=it.next();
			double[] dis=(double[])h1.get(o);
			h3.put(o, dis);
		}
		
		s=h2.keySet();
	    it=s.iterator();
		while(it.hasNext()){
			Object o=it.next();
			double[] dis=(double[])h2.get(o);
			h3.put(o, dis);
		}
		
		return h3;
	}
	
	//In use 
	public ArrayList getInnernodesNameList(HashMap MDS_Name_Dis){
		ArrayList al=new ArrayList();
		
		Set s=MDS_Name_Dis.keySet();
		Iterator it=s.iterator();
		while(it.hasNext()){
			String name=(String)it.next();
			al.add(name);
		}
		
		return al;
	}
	
	
	//In use 
	//generate leaf name coordinates 
	public void textDigit(HashMap Leaf_Name_Coords,database db, File fout){
			
		Visual vis=new Visual();
		StringBuffer sb=new StringBuffer();
		StringBuffer sb2=new StringBuffer();
		Set s=Leaf_Name_Coords.keySet();
		Iterator it=s.iterator();
		int cnt=0;
		while(it.hasNext()){
			it.next();
			cnt++;
		}
		int coordsize=cnt;
		double[][] coords=new double[cnt][2];
		ArrayList lstrainname=new ArrayList();
		it=s.iterator();
		cnt=1;
		while(it.hasNext()){
			String ss=(String)it.next();
			double[] coord=(double[])Leaf_Name_Coords.get(ss);
			double[] newcoord=vis.nonOverlapcoords(coords, coord, cnt-1);
			coords[cnt-1][0]=newcoord[0];
			coords[cnt-1][1]=newcoord[1];
			String strainname=db.findByID(ss).getAnnotation();
			lstrainname.add(strainname);
			//sb.append(cnt+"	"+(newcoord[0]+0.001)+"	"+(newcoord[1]+0.001)+"\n");
			//sb2.append(cnt+":"+strainname.substring(1,s_l-1)+"\n");
			cnt++;
		}
		
		boolean ischange=true;
		while(ischange){
			ischange=false;
			for(int i=0;i<(coordsize-1);i++){
				double a=coords[i][0];
				double b=coords[i+1][0];
				if(a>b){
					double[] coordc=coords[i];
					coords[i]=coords[i+1];
					coords[i+1]=coordc;
					
					String sa=(String)lstrainname.get(i);
					String sB=(String)lstrainname.get(i+1);
					String sc=sa;
					lstrainname.set(i, sB);
					lstrainname.set(i+1,sc);
					ischange=true;
				}
			}
		}
		
		
		for(int i=0;i<coordsize;i++){
			sb.append((i+1)+"	"+(coords[i][0]+0.001)+"	"+(coords[i][1]+0.001)+"\n");
			sb2.append((i+1)+":"+(String)lstrainname.get(i)+"\n");
		}
		
		
		try{
			FileWriter fw=new FileWriter(fout);
			fw.write(sb.toString());
			fw.close();
			
			File f2=new File(fout.getAbsolutePath()+".content");
			fw=new FileWriter(f2);
			fw.write(sb2.toString());
			fw.close();
		}catch(Exception e){
			e.printStackTrace();
		}
	}
	
	//In use 
	//generate leaf name coordinates and mark the selected sequences from eplist
	public void textDigit(HashMap Leaf_Name_Coords, File fout, List eplist, database db){

		Visual vis=new Visual();
		StringBuffer sb=new StringBuffer();
		StringBuffer sb2=new StringBuffer();
		StringBuffer sb3=new StringBuffer();
		Set s=Leaf_Name_Coords.keySet();
		Iterator it=s.iterator();
		int cnt=0;
		while(it.hasNext()){
			it.next();
			cnt++;
		}
		int coordsize=cnt;
		double[][] coords=new double[cnt][2];
		ArrayList lstrainname=new ArrayList();
		ArrayList lacc=new ArrayList();
		it=s.iterator();
		cnt=1;
		while(it.hasNext()){
			String ss=(String)it.next();
			lacc.add(ss);
			double[] coord=(double[])Leaf_Name_Coords.get(ss);
			double[] newcoord=vis.nonOverlapcoords(coords, coord, cnt-1);
			coords[cnt-1][0]=newcoord[0];
			coords[cnt-1][1]=newcoord[1];
			String strainname=db.findByID(ss).getAnnotation();
			lstrainname.add(strainname);
			
			//sb.append(cnt+"	"+(newcoord[0]+0.001)+"	"+(newcoord[1]+0.001)+"\n");
			//sb2.append(cnt+":"+strainname.substring(1,s_l-1)+"\n");
			cnt++;
		}
		
		boolean ischange=true;
		while(ischange){
			ischange=false;
			for(int i=0;i<(coordsize-1);i++){
				double a=coords[i][0];
				double b=coords[i+1][0];
				if(a>b){
					double[] coordc=coords[i];
					coords[i]=coords[i+1];
					coords[i+1]=coordc;
					
					String sa=(String)lstrainname.get(i);
					String sB=(String)lstrainname.get(i+1);
					String sc=sa;
					lstrainname.set(i, sB);
					lstrainname.set(i+1,sc);
					
					sa=(String)lacc.get(i);
					sB=(String)lacc.get(i+1);
					sc=sa;
					lacc.set(i, sB);
					lacc.set(i+1, sc);
					
					ischange=true;
				}
			}
		}
		
		
		for(int i=0;i<coordsize;i++){
			String acc=(String)lacc.get(i);
			if(eplist.contains(acc)){
				sb.append((i+1)+"	"+(coords[i][0]+0.0005)+"	"+(coords[i][1]+0.0005)+"	"+"1"+"\n");
				sb2.append((i+1)+":"+(String)lstrainname.get(i)+" *\n");
				sb3.append(acc+"	"+(i+1)+". "+(String)lstrainname.get(i)+" *\n");
			}else{
				sb.append((i+1)+"	"+(coords[i][0]+0.001)+"	"+(coords[i][1]+0.001)+"	"+"0"+"\n");
				sb2.append((i+1)+":"+(String)lstrainname.get(i)+"\n");
				sb3.append(acc+"	"+(i+1)+". "+(String)lstrainname.get(i)+"\n");
			}
		}
		
		
		try{
			FileWriter fw=new FileWriter(fout);
			fw.write(sb.toString());
			fw.close();
			
			File f2=new File(fout.getAbsolutePath()+".content");
			fw=new FileWriter(f2);
			fw.write(sb2.toString());
			fw.close();
			
			File f3=new File(fout.getAbsolutePath()+".help");
			fw=new FileWriter(f3);
			fw.write(sb3.toString());
			fw.close();
		}catch(Exception e){
			e.printStackTrace();
		}
	}
	 
	
	//In use 
	public HashMap runPnp(File fdb, File fTextout, File treefile, File fLeafCoords, int numPoint, double xlim1, double xlim2, double ylim1, double ylim2, int maxiters, File fcc, File fcc_bezier,double scale_factor) throws Exception{
		
		numPoint=numPoint-2;
		
		//Parse the tree and get tree root
		NewickParse parser = new NewickParse(); 
		node root=parser.parseFromFile(treefile);
		this.tree_root=root;
		
		//Correct minus distances
		postOrder po=new postOrder();
		po.correctMinusDisInTree(root);
		
		//fill in node_list and tree distance matrix and calculating C
		po.setScaleFactor(scale_factor);
		double C=po.sumTreeSpacedistance(root);
		
		//load Leaf name to coords
		HashMap Leaf_Name_Coords=this.loadLeafCoords(fLeafCoords);
		
		//this.text(Leaf_Name_Coords, fPair, fTextout);
		//generate leaf name coordinates 
		database db=new database(fdb);
		//this.textDigit(Leaf_Name_Coords,db, fTextout);
		this.textDigit(Leaf_Name_Coords, fTextout, new LinkedList(), db);
		
		
		//Inner nodes name to 3 connected nodes distances
		po.postTree_cal_tree_dis(root.getLeft());
		HashMap Tree_Name_Dis=po.getTreeNameToDis();
		
		//Above calculated from the tree
		////////////////////////////////////////////////////////////////////////////
		//Below calculate the inner nodes
		
		//Init inner nodes coords and assign to inner nodes names
		double[][] coords=this.init(numPoint, xlim1, xlim2, ylim1, ylim2);
		HashMap Tree_name_Coords=this.assignCoordsToInnernodes(Tree_Name_Dis, coords);
		
		//Vertices name to coords, name and double[] pair
		HashMap Name_Coords=this.intergateHMs(Leaf_Name_Coords, Tree_name_Coords);
		
		//fill in mds distance matrix, this matrix will be used for shannon mapping 
		po.cal_mds_dis(root, Name_Coords);
		
		//Cal intit error 
		double error0=po.calError(C);
		double error=error0;
		System.out.println("init error:"+error0);
		
		//To store the best solution find during calculation
		double smallError=Double.POSITIVE_INFINITY;
		HashMap small_Name_Coords=new HashMap();
		
		//Inner nodes name list
		ArrayList InnerNodesNameList=this.getInnernodesNameList(Tree_Name_Dis);
		
		//All nodes: inner and leaf nodes list
		ArrayList NodesNameList=po.getNodeList();
		
		
		ArrayList al=null;
		double cc=0;

		for(int i=0;i<maxiters;i++){
			
			int[] randidx=Permutation.next(numPoint);
			for(int k=0;k<numPoint;k++){
				//Select one innder node and update it against all nodes including leaf nodes 
				int idx=randidx[k];
				String tname=(String)InnerNodesNameList.get(idx);
				
				//Update the selected point
				//Execute the normal update for all nodes once and update the nodes with edges smaller than the in the tree 4 times
				if(i % 5==0){
					//Name_Coords=po.upDateSmall(tname, Name_Coords);
					//Name_Coords=po.upDateBelow(tname, Name_Coords);
					Name_Coords=po.upDate(tname, Name_Coords);
					
				}else{
					Name_Coords=po.upDateBelow(tname, Name_Coords);
					//Name_Coords=po.upDateSmall(tname, Name_Coords);
				}
				
				
				//Name_Coords=po.upDate(tname, Name_Coords);
				//Name_Coords=po.upDateSmall(tname, Name_Coords);
				
				//Calculate new MDS distances
				po.cal_mds_dis(root, Name_Coords);
				
			}
			al=po.ccMatrix();
			cc=po.Calcc(al);
			error=po.calError(C);
			
			if (error<smallError){
				smallError=error;
				small_Name_Coords=(HashMap)Name_Coords.clone();
			}
			
			//System.out.println("Error in Iters "+i+": "+error);
		}
		
		//System.out.println("init error:"+error0);
		//small_Name_Coords=Name_Coords;
		po.cal_mds_dis(root, small_Name_Coords);
		po.ccPlot(fcc);
		al=po.ccMatrix();
		cc=po.Calcc(al);
		
		
		po.calBezier(root, small_Name_Coords, 1);
		this.hm_controlpoint=po.getControlp();
		po.ccPlotBezier(fcc_bezier);
		ArrayList al_bezier=po.ccMatrixBezier();
		double cc_bezier=po.Calcc(al_bezier);
		double error_bezier=po.calErrorBezier(C);

		System.out.println("sacle:"+scale_factor+"	CC:"+cc+"	Bezier CC:"+cc_bezier+"	Error:"+smallError+"	Error after Bezier:"+error_bezier);
		
		File ferror=new File(fcc_bezier.getAbsolutePath()+".error");
		
		try{
			FileWriter fw=new FileWriter(ferror);
			fw.write("sacle:"+scale_factor+"	CC:"+cc+"	Bezier CC:"+cc_bezier+"	Error:"+smallError+"	Error after Bezier:"+error_bezier);
			fw.close();
		}catch(Exception e){
			e.printStackTrace();
		}
		
		this.sammon_error=smallError;
		this.bezier_error=error_bezier;
		return small_Name_Coords;
	}
	
	
	//In use 
	public HashMap runPnp(database db, File fTextout, File treefile, File fLeafCoords, int numPoint, double xlim1, double xlim2, double ylim1, double ylim2, int maxiters, File fcc, File fcc_bezier,double scale_factor) throws Exception{
		
		numPoint=numPoint-2;
		
		//Parse the tree and get tree root
		NewickParse parser = new NewickParse(); 
		node root=parser.parseFromFile(treefile);
		this.tree_root=root;
		
		//Correct minus distances
		postOrder po=new postOrder();
		po.correctMinusDisInTree(root);
		
		//fill in node_list and tree distance matrix and calculating C
		po.setScaleFactor(scale_factor);
		double C=po.sumTreeSpacedistance(root);
		
		//load Leaf name to coords
		HashMap Leaf_Name_Coords=this.loadLeafCoords(fLeafCoords);
		
		//this.text(Leaf_Name_Coords, fPair, fTextout);
		//generate leaf name coordinates 
		//database db=new database(fdb);
		//this.textDigit(Leaf_Name_Coords,db, fTextout);
		this.textDigit(Leaf_Name_Coords, fTextout, new LinkedList(), db);
		
		
		//Inner nodes name to 3 connected nodes distances
		po.postTree_cal_tree_dis(root.getLeft());
		HashMap Tree_Name_Dis=po.getTreeNameToDis();
		
		//Above calculated from the tree
		////////////////////////////////////////////////////////////////////////////
		//Below calculate the inner nodes
		
		//Init inner nodes coords and assign to inner nodes names
		double[][] coords=this.init(numPoint, xlim1, xlim2, ylim1, ylim2);
		HashMap Tree_name_Coords=this.assignCoordsToInnernodes(Tree_Name_Dis, coords);
		
		//Vertices name to coords, name and double[] pair
		HashMap Name_Coords=this.intergateHMs(Leaf_Name_Coords, Tree_name_Coords);
		
		//fill in mds distance matrix, this matrix will be used for shannon mapping 
		po.cal_mds_dis(root, Name_Coords);
		
		//Cal intit error 
		double error0=po.calError(C);
		double error=error0;
		System.out.println("init error:"+error0);
		
		//To store the best solution find during calculation
		double smallError=Double.POSITIVE_INFINITY;
		HashMap small_Name_Coords=new HashMap();
		
		//Inner nodes name list
		ArrayList InnerNodesNameList=this.getInnernodesNameList(Tree_Name_Dis);
		
		//All nodes: inner and leaf nodes list
		ArrayList NodesNameList=po.getNodeList();
		
		
		ArrayList al=null;
		double cc=0;
		
		//todo: potential bugs
		numPoint = InnerNodesNameList.size();
		for(int i=0;i<maxiters;i++){
			
			int[] randidx=Permutation.next(numPoint);
			for(int k=0;k<numPoint;k++){
				//Select one innder node and update it against all nodes including leaf nodes 
				int idx=randidx[k];
				String tname=(String)InnerNodesNameList.get(idx);
				
				//Update the selected point
				//Execute the normal update for all nodes once and update the nodes with edges smaller than the in the tree 4 times
				if(i % 5==0){
					//Name_Coords=po.upDateSmall(tname, Name_Coords);
					//Name_Coords=po.upDateBelow(tname, Name_Coords);
					Name_Coords=po.upDate(tname, Name_Coords);
					
				}else{
					Name_Coords=po.upDateBelow(tname, Name_Coords);
					//Name_Coords=po.upDateSmall(tname, Name_Coords);
				}
				
				
				//Name_Coords=po.upDate(tname, Name_Coords);
				//Name_Coords=po.upDateSmall(tname, Name_Coords);
				
				//Calculate new MDS distances
				po.cal_mds_dis(root, Name_Coords);
				
			}
			al=po.ccMatrix();
			cc=po.Calcc(al);
			error=po.calError(C);
			
			if (error<smallError){
				smallError=error;
				small_Name_Coords=(HashMap)Name_Coords.clone();
			}
			
			//System.out.println("Error in Iters "+i+": "+error);
		}
		
	
		//small_Name_Coords=Name_Coords;
		po.cal_mds_dis(root, small_Name_Coords);
		po.ccPlot(fcc);
		al=po.ccMatrix();
		cc=po.Calcc(al);
		
		po.calBezier(root, small_Name_Coords, 1);
		this.hm_controlpoint=po.getControlp();
		po.ccPlotBezier(fcc_bezier);
		ArrayList al_bezier=po.ccMatrixBezier();
		double cc_bezier=po.Calcc(al_bezier);
		double error_bezier=po.calErrorBezier(C);

		System.out.println("Mapping results:");
		System.out.println("scale factor:	"+scale_factor);
		System.out.println("correlation coefficient:	"+cc);
		System.out.println("correlation coefficient after bezier compensation:	"+cc_bezier);
		System.out.println("error:	"+ smallError);
		System.out.println("final error after bezier compensation:	"+error_bezier);
		
		File ferror=new File(fcc_bezier.getAbsolutePath()+".error");
		
		try{
			FileWriter fw=new FileWriter(ferror);
			fw.write("sacle:"+scale_factor+"	CC:"+cc+"	Bezier CC:"+cc_bezier+"	Error:"+smallError+"	Error after Bezier:"+error_bezier);
			fw.close();
		}catch(Exception e){
			e.printStackTrace();
		}
		
		this.sammon_error=smallError;
		this.bezier_error=error_bezier;
		return small_Name_Coords;
	}

	
	//In use. Marking the selcted sequences using the eplist 
	public HashMap runPnp(List eplist, File fdb, File fTextout, File treefile, File fLeafCoords, int numPoint, double xlim1, double xlim2, double ylim1, double ylim2, int maxiters, File fcc, File fcc_bezier,double scale_factor) throws Exception{
		
		numPoint=numPoint-2;
		
		//Parse the tree and get tree root
		NewickParse parser = new NewickParse(); 
		node root=parser.parseFromFile(treefile);
		this.tree_root=root;
		
		//Correct minus distances
		postOrder po=new postOrder();
		po.correctMinusDisInTree(root);
		
		//fill in node_list and tree distance matrix and calculating C
		po.setScaleFactor(scale_factor);
		double C=po.sumTreeSpacedistance(root);
		
		//load Leaf name to coords
		HashMap Leaf_Name_Coords=this.loadLeafCoords(fLeafCoords);
		
		//generate leaf name coordinates
		database db=new database(fdb);
		this.textDigit(Leaf_Name_Coords, fTextout, eplist, db);
		
		//Inner nodes name to 3 connected nodes distances
		po.postTree_cal_tree_dis(root.getLeft());
		HashMap Tree_Name_Dis=po.getTreeNameToDis();
		
		//Above calculated from the tree
		////////////////////////////////////////////////////////////////////////////
		//Below calculate the inner nodes
		
		//Init inner nodes coords and assign to inner nodes names
		double[][] coords=this.init(numPoint, xlim1, xlim2, ylim1, ylim2);
		HashMap Tree_name_Coords=this.assignCoordsToInnernodes(Tree_Name_Dis, coords);
		
		//Vertices name to coords, name and double[] pair
		HashMap Name_Coords=this.intergateHMs(Leaf_Name_Coords, Tree_name_Coords);
		
		//fill in mds distance matrix, this matrix will be used for shannon mapping 
		po.cal_mds_dis(root, Name_Coords);
		
		//Cal intit error 
		double error0=po.calError(C);
		double error=error0;
		System.out.println("init error:"+error0);
		
		//To store the best solution find during calculation
		double smallError=Double.POSITIVE_INFINITY;
		HashMap small_Name_Coords=new HashMap();
		
		//Inner nodes name list
		ArrayList InnerNodesNameList=this.getInnernodesNameList(Tree_Name_Dis);
		
		//All nodes: inner and leaf nodes list
		ArrayList NodesNameList=po.getNodeList();
		
		ArrayList al=null;
		double cc=0;

		for(int i=0;i<maxiters;i++){
			
			int[] randidx=Permutation.next(numPoint);
			for(int k=0;k<numPoint;k++){
				//Select one innder node and update it against all nodes including leaf nodes 
				int idx=randidx[k];
				String tname=(String)InnerNodesNameList.get(idx);
				
				//Update the selected point
				//Execute the normal update for all nodes once and update the nodes with edges smaller than the in the tree 4 times
				if(i % 5==0){
					Name_Coords=po.upDate(tname, Name_Coords);
				}else{
					Name_Coords=po.upDateBelow(tname, Name_Coords);
				}
				
				//Calculate new MDS distances
				po.cal_mds_dis(root, Name_Coords);
				
			}
			al=po.ccMatrix();
			cc=po.Calcc(al);
			error=po.calError(C);
			
			if (error<smallError){
				smallError=error;
				small_Name_Coords=(HashMap)Name_Coords.clone();
			}
			
			System.out.println("C:"+C+","+"error in Iters "+i+": "+error);
		}
		
		po.cal_mds_dis(root, small_Name_Coords);
		po.ccPlot(fcc);
		al=po.ccMatrix();
		cc=po.Calcc(al);
		
		
		po.calBezier(root, small_Name_Coords, 1);
		this.hm_controlpoint=po.getControlp();
		po.ccPlotBezier(fcc_bezier);
		ArrayList al_bezier=po.ccMatrixBezier();
		double cc_bezier=po.Calcc(al_bezier);
		double error_bezier=po.calErrorBezier(C);

		System.out.println("sacle:"+scale_factor+"	CC:"+cc+"	Bezier CC:"+cc_bezier+"	Error:"+smallError+"	Error after Bezier:"+error_bezier);
		
		File ferror=new File(fcc_bezier.getAbsolutePath()+".error");
		
		try{
			FileWriter fw=new FileWriter(ferror);
			fw.write("sacle:"+scale_factor+"	CC:"+cc+"	Bezier CC:"+cc_bezier+"	Error:"+smallError+"	Error after Bezier:"+error_bezier);
			fw.close();
		}catch(Exception e){
			e.printStackTrace();
		}
		
		this.sammon_error=smallError;
		this.bezier_error=error_bezier;
		return small_Name_Coords;
	}
	
	public HashMap copyHm(HashMap hm){
		HashMap new_hm=new HashMap();
		Set s=hm.keySet();
		Iterator it=s.iterator();
		while(it.hasNext()){
			Object name=it.next();
			Object o=hm.get(name);
			new_hm.put(name, o);
		}
		return new_hm;
	}
	
	public double maxScale(HashMap hmTreeDis, HashMap hmMDSDis){
		double t=0;
		Set s=hmTreeDis.keySet();
		Iterator it=s.iterator();
		while(it.hasNext()){
			String name=(String)it.next();
			Double Dtree=(Double)hmTreeDis.get(name);
			Double Dmds=(Double)hmMDSDis.get(name);
			double t_temp=Dmds.doubleValue()/Dtree.doubleValue();
			if(t_temp>t){
				t=t_temp;
			}
		}
		
		return t;
	}
	
	public double LBezier(double x0, double y0, double x1, double y1, double x2, double y2){
		double length=0;
		
		double ax=x0-2*x1+x2;
		double ay=y0-2*y1+y2;
		double bx=2*x1-2*x0;
		double by=2*y1-2*y0;
		double A=4*(ax*ax+ay*ay);
		double B=4*(ax*bx+ay*by);
		double C=bx*bx+by*by;
		double A32=Math.pow(A, 3.0/2.0);
		double SABC=Math.sqrt(A+B+C);
		double SA=Math.sqrt(A);
		double SC=Math.sqrt(C);
		
		length=(4*A32*SABC+2*SA*B*(SABC-SC)+(4*C*A-B*B)*Math.log((2*SA+B/SA+2*SABC)/((B/SA+2*SC))))/(8*A32);
		
		return length;
		
		
	}
	
	public double dis(double[] c1, double[] c2){
		double dis=0;
		for(int i=0;i<c1.length;i++){
			dis+=((c1[i]-c2[i])*(c1[i]-c2[i]));
		}
		return Math.sqrt(dis);
	}
	
	public double[] findControlPoint(double x0, double y0, double x2, double y2, double scale, double tree_dis, double error){
		double[] control =new double[3];
		double[] c1=new double[2];
		c1[0]=x0;
		c1[1]=y0;
		double[] c2=new double[2];
		c2[0]=x2;
		c2[1]=y2;
		double Lcord=this.dis(c1, c2);
		double Ltarget=scale*tree_dis;
		double Lstring=Lcord;
		
		//double k=-(y2-y0)/(x2-x0);
		double k=-(x2-x0)/(y2-y0);
		double xm=(x0+x2)/2.0;
		double ym=(y0+y2)/2.0;
		
		double step=0.00001;
		double x1=xm;
		double y1=k*(x1-xm)+ym;
		
		if(Lcord>=Ltarget){
			control[0]=x0;
			control[1]=y0;
			control[2]=Lcord;
		}else{
		
			double c_error=Ltarget-Lstring;
			while(c_error>0){
				x1=x1+step;
				y1=k*(x1-xm)+ym;
				Lstring=this.LBezier(x0, y0, x1, y1, x2, y2);
				c_error=Ltarget-Lstring;
			}
			double abs_c_error=Math.sqrt(c_error*c_error);
			while(abs_c_error>error){
				step=step/2;
				if(c_error>0){
					x1=x1+step;
					y1=k*(x1-xm)+ym;
					Lstring=this.LBezier(x0, y0, x1, y1, x2, y2);
					c_error=Ltarget-Lstring;
					
				}else{
					x1=x1-step;
					y1=k*(x1-xm)+ym;
					Lstring=this.LBezier(x0, y0, x1, y1, x2, y2);
					c_error=Ltarget-Lstring;
				}
				abs_c_error=Math.sqrt(c_error*c_error);
			}
			
			control[0]=x1;
			control[1]=y1;
			control[2]=Lstring;
		}
		
		//System.out.println("oL:"+Lcord+"	tL:"+Ltarget+"	bL:"+Lstring);
		return control;
		
	}
	
	//In use 
	//Generate Innernodes coordinates, will be plot to distinguish the leaf
	public void Name_Coords_toString(HashMap name_coord, File fncout) throws Exception{
		StringBuffer sb=new StringBuffer();
		Set s=name_coord.keySet();
		Iterator it=s.iterator();
		while(it.hasNext()){
			String o=(String)it.next();
			if(o.startsWith("#")){
				double[] coord=(double[])name_coord.get(o);
				sb.append(coord[0]+"	"+coord[1]+"\n");
			}
		}
		FileWriter fw=new FileWriter(fncout);
		fw.write(sb.toString());
		fw.close();
	}
	
	public static void main(String[] args) {
		


	}

}

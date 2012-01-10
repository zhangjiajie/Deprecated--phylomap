package pal.tree;
import java.io.*;
import java.util.*;
import pal.misc.*;

public class ExtractSubtree {

	/**
	 * @param args
	 */
	private SimpleTree tFull;
	private SimpleTree tSub;
	private HashSet hsSubtreeNodes=new HashSet();
	private double pd=0;
	
	public ExtractSubtree (SimpleTree tFull, HashSet hsSubtreeNodes){
		this.tFull=tFull;
		this.hsSubtreeNodes=hsSubtreeNodes;
		//it is very important to id the nodes as internal or leaf in the beginning
		for(int i=0;i<tFull.getExternalNodeCount();i++){
			SimpleNode enote=(SimpleNode)tFull.getExternalNode(i);
			enote.setIsLeaf(true);
		}
		
		for(int i=0;i<tFull.getInternalNodeCount();i++){
			SimpleNode inote=(SimpleNode)tFull.getInternalNode(i);
			inote.setIsLeaf(false);
			Identifier id= inote.getIdentifier();
			id.setName("inode"+i);
			inote.setIdentifier(id);
		}
		
		//System.out.println("Read in the tree: \n"+tFull.toString());
		
	}
	
	public ExtractSubtree (String sFulltreeFilename, HashSet hsSubtreeNodes){
		try{
			tFull= new ReadTree(sFulltreeFilename);
			this.hsSubtreeNodes=hsSubtreeNodes;
			
			//it is very important to id the nodes as internal or leaf in the beginning
			for(int i=0;i<tFull.getExternalNodeCount();i++){
				SimpleNode enote=(SimpleNode)tFull.getExternalNode(i);
				enote.setIsLeaf(true);
			}
			
			for(int i=0;i<tFull.getInternalNodeCount();i++){
				SimpleNode inote=(SimpleNode)tFull.getInternalNode(i);
				inote.setIsLeaf(false);
				Identifier id= inote.getIdentifier();
				id.setName("inode"+i);
				inote.setIdentifier(id);
			}
			
			//System.out.println("Read in the tree: \n"+tFull.toString());
		}catch(Exception e){
			e.printStackTrace();
		}
	}
	
	public SimpleTree getSubtree(){
		//System.out.println(tFull.toString());
		if(!this.isBinary(tFull)){
			System.out.println("The input tree is not binary!!");
		}
		//if(this.isBinary(tFull)){
		if(true){	
			//find the leaf to delete 
			LinkedList ldel=new LinkedList();
			for(int i=0; i<tFull.getExternalNodeCount();i++){
				Node node=tFull.getExternalNode(i);
				if(!hsSubtreeNodes.contains(node.getIdentifier().getName())){
					ldel.add(node.getIdentifier().getName());
				}
			}
			
			//delete one leaf at a time
			for(int i=0; i<ldel.size(); i++){
				String nname=(String)ldel.get(i);
				this.postOrderdel(tFull.getRoot(), nname);
				//this.cleanTree(tFull);
			}
			//this.lastOrder(tFull.getRoot());
			this.tSub=(SimpleTree)tFull.getCopy();
			//System.out.println("tree before clean: \n"+tSub.toString());
			
			//cleaning the subtree
			this.cleanTree(this.tFull);
			//System.out.println("tree after clean: \n"+tFull.toString());
			
			//only run this for the final step!!
			//this.removeInnerNodename(tSub);
			this.removeInnerNodename(tFull);
			//System.out.println("messy tree after remove innernode name: \n"+tSub.toString());
			//System.out.println("clean tree after remove innernode name: \n"+tFull.toString());
			return this.tFull;
		}else{
			System.out.println("The input tree is not binary!!");
			return null;
		}
	}
	
	public static void writeSubTree(File fout, SimpleTree st){
		try{
			FileWriter fw=new FileWriter(fout);
			fw.write(st.toString());
			fw.close();
		}catch(Exception e){
			e.printStackTrace();
		}
	}
	
	//bugs, do not use!
	//delete the nodes that are not in the sub notes list
	private void lastOrder(Node sn){
		if(sn.isLeaf()){
			if(!hsSubtreeNodes.contains(sn.getIdentifier().getName())){
				Node sp=sn.getParent();
				for(int i=0; i<sp.getChildCount(); i++){
					if(!hsSubtreeNodes.contains(sp.getChild(i).getIdentifier().getName())){
						sp.removeChild(i);
						//System.out.println("renoved:" + sdel);
					}
				}
			}
			
		}else{
			for(int i=0; i<sn.getChildCount();i++){
				this.lastOrder(sn.getChild(i));
			}
		}
	}
	
	//delete one leaf at a time
	private void postOrderdel(Node sn, String nodename){
		
		for(int i=0; i<sn.getChildCount();i++){
			this.postOrderdel(sn.getChild(i), nodename);
		}
		
		SimpleNode ssn=(SimpleNode)sn;
		
		if(ssn.getIsLeaf()){
			if(nodename.equals(sn.getIdentifier().getName())){
				Node sp=sn.getParent();
				if(sp==null){
					System.out.println("Warnning:: parent is null!! (in postOrderdel)");
				}
				for(int j=0; j<sp.getChildCount(); j++){
					if(nodename.equals(sp.getChild(j).getIdentifier().getName())){
						Node rnode=sp.removeChild(j);
						rnode.setParent(null);
						//System.out.println("removed:"+ rnode.getIdentifier().getName());
						break;
					}
				}
			}
		}
	}
	
	private boolean hasChanged=false;
	//postOrder clean the tree, will only do one operation for each call 
	private void treeStructureClean(Node root){
		for(int i=0;i<root.getChildCount();i++){
			this.treeStructureClean(root.getChild(i));
		}
		SimpleNode sn= (SimpleNode)root;
		if(sn.getIsLeaf()){
			//do nothing!
			//System.out.println("Warning:: reached a leaf node: "+root.getIdentifier().getName());
		}else{
			if(root.isRoot()){
				if(root.getChildCount()==1 && !sn.getIsLeaf()){
					//System.out.println(root.getIdentifier().getName());
					this.tFull.setRoot(root.getChild(0));
					root=root.getChild(0);
					root.setParent(null);
					hasChanged=true;
					//System.out.println(root.getIdentifier().getName());
					//System.out.println("Operation:: move root down");
				}
			}else if(root.getChildCount()==0 && !hasChanged){
				Node parent = root.getParent();
				if(parent!=null){
					for(int i=0;i<parent.getChildCount();i++){
						if(root.getIdentifier().getName().equals(parent.getChild(i).getIdentifier().getName())){
							parent.removeChild(i);
							break;
						}
					}
					hasChanged=true;
					//System.out.println("Operation:: del one inner node");
				}
			}else if(root.getChildCount()==1 && !hasChanged){
				Node parent = root.getParent();
				Node child = root.getChild(0);
				if(parent!=null){
					child.setBranchLength(child.getBranchLength()+root.getBranchLength());
					for(int i=0;i<parent.getChildCount();i++){
						if(root.getIdentifier().getName().equals(parent.getChild(i).getIdentifier().getName())){
							parent.removeChild(i);
							parent.addChild(child);
							break;
						}
					}
					hasChanged=true;
					//System.out.println("Operation:: move one inner node");
				}
			}
		}
	}
	
	//test if the tree is binary tree 
	public boolean isBinary(SimpleTree st){
		boolean flag=true;
		for(int i=0;i<st.getInternalNodeCount();i++){
			Node inode=st.getInternalNode(i);
			if(inode.getChildCount()!=2 && !inode.isRoot()){
				flag=false;
				//System.out.println (inode.getChildCount());
				//System.out.println (inode.getIdentifier().getName());
				//System.out.println (inode.getChild(0).getIdentifier().getName());
				//System.out.println (inode.getChild(1).getIdentifier().getName());
				//System.out.println (inode.getChild(2).getIdentifier().getName());
				//System.out.println ();
			}
		}
		return flag;
	}
	
	//clean the tree 
	private void cleanTree(SimpleTree st){
			this.hasChanged=true;
			while(this.hasChanged){
				this.hasChanged=false;
				this.treeStructureClean(st.getRoot());
				//System.out.println(st.toString());
			}
	}
	
	private void removeInnerNodename(SimpleTree st){
		st.createNodeList();
		for(int i=0;i<st.getInternalNodeCount();i++){
			Node inode=st.getInternalNode(i);
			Identifier id=inode.getIdentifier();
			id.setName("");
			inode.setIdentifier(id);
		}
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
	
	public static void main(String[] args) {
		// TODO Auto-generated method stub
		String treefile="/home/zhangje/data/Tree/outtree";
		File fout=new File("/home/zhangje/data/Tree/subtree.tre");
		HashSet hs=new HashSet();
		hs.add("ABD77822");
		hs.add("ABS00325");
		hs.add("ACF54514");
		hs.add("ACA24583");
		hs.add("AAG01204");
		hs.add("AAV48837");
		hs.add("AAA43482");
		hs.add("ABW86435");
		hs.add("ACD35866");
		ExtractSubtree es=new ExtractSubtree(treefile, hs);
		SimpleTree st=es.getSubtree();
	    if (st!=null){
	    		es.writeSubTree(fout, st);
	    }else{
	    		System.out.println("tree is null");
	    }
	}

}

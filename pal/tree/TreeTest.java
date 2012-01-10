package pal.tree;
import java.io.*;

public class TreeTest {

	/**
	 * @param args
	 */
	public void lastOrder(Node sn, String sdel){
		if(sn.isLeaf()){
			SimpleNode ssn=(SimpleNode)sn;
			System.out.println(ssn.getIdentifier().getName()+":"+ssn.getBranchLength()+":"+ssn.getBranchLengthSE());
			if(sdel.equals(ssn.getIdentifier().getName())){
				SimpleNode sp=(SimpleNode)ssn.getParent();
				for(int i=0; i<sp.getChildCount(); i++){
					if(sdel.equals(sp.getChild(i).getIdentifier().getName())){
						sp.removeChild(i);
						System.out.println("renoved:" + sdel);
					}
				}
			}
			
		}else{
			if(sn.isRoot()){
				System.out.println("root"+sn.getIdentifier().getName()+":"+sn.getBranchLength()+":"+sn.getBranchLengthSE());
				System.out.println ("rootparent:"+sn.getParent());
			}
			for(int i=0; i<sn.getChildCount();i++){
				this.lastOrder(sn.getChild(i), sdel);
			}
		}
	}
	
	
	public void delNote(String name, String sftree){
		try{
			ReadTree tree1=new ReadTree(sftree);
			System.out.println(tree1.toString());
			this.lastOrder(tree1.getRoot(), name);
			System.out.println(tree1.toString());
			ReadTree tree2=new ReadTree(sftree);
			
			TreeUtils tu=new TreeUtils();
			//double drf=tu.getRobinsonFouldsDistance(tree1, tree2);
			//System.out.println(drf);
		}catch(Exception e){
			e.printStackTrace();
		}
		
	}
	
	public static void main(String[] args) {
		// TODO Auto-generated method stub
		String sfile="/home/zhangje/data/Tree/outtree";
		TreeTest tt=new TreeTest();
		tt.delNote("ABD77822", sfile);
	}

}

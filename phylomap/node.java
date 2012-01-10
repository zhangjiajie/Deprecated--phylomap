package phylomap;

import java.io.Serializable;

public class node implements Serializable{

	/**
	 * @param args
	 */
	
	double distance=0;
	node left=null;
	node right=null;
	String id=null;
	int nodecnt=1;
	String info=null;
	
	
	public void setInfo(String info){
		this.info=info;
	}
	
	public String getInfo()
	{
		return this.info;
	}
	public void setLeft(node left){
		this.left=left;
	}
	public node getLeft(){
		return this.left;
	}
	public void setRight(node right){
		this.right=right;
	}
	public node getRight(){
		return this.right;
	}
	public void setId(String id){
		this.id=id;
	}
	public String getId(){
		return this.id;
	}
	
	public double getDistance(){
		return this.distance;
	}
	
	public void setDistance(double distance){
		this.distance=distance;
	}
	
	public int getNodecnt(){
		return this.nodecnt;
	}
	
	public void setNodecnt(){
		if(left!=null && right!=null){
			this.nodecnt=left.getNodecnt()+right.getNodecnt();
		}
	}
	
	
	
	
	public static void main(String[] args) {
		// TODO Auto-generated method stub

	}

}



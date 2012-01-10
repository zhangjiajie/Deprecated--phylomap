package pal.sampling;
import pal.tree.*;


public class TreeDistance {

	/**
	 * @param args
	 */
	public static double RFdistance(String t1, String t2){
		double dis=-1;
		try{
			ReadTree tree1 = new ReadTree(t1);
			ReadTree tree2 = new ReadTree(t2);
			dis=TreeUtils.getRobinsonFouldsDistance(tree1, tree2);
		}catch (Exception e){
			e.printStackTrace();
		}
		return dis;
	}
	
	public static double RFdistanceScaled(String t1, String t2){
		double dis=-1;
		try{
			ReadTree tree1 = new ReadTree(t1);
			ReadTree tree2 = new ReadTree(t2);
			dis=TreeUtils.getRobinsonFouldsRescaledDistance(tree1, tree2);
		}catch (Exception e){
			//e.printStackTrace();
		}
		return dis;
	}
	
	public static double RFdistanceScaled(ReadTree t1, ReadTree t2){
		double dis=-1;
		try{
			dis=TreeUtils.getRobinsonFouldsRescaledDistance(t1, t2);
		}catch (Exception e){
			//e.printStackTrace();
		}
		return dis;
	}
	
	
	public static void main(String[] args) {
		// TODO Auto-generated method stub
		//String t1="/home/zhangje/data/alexi/d4114/randsample/rand102.tre";
		//String t2="/home/zhangje/storage/data/d4114/randsample/randsample/RAxML_bestTree.rand102";
		
		String t1="/home/zhangje/storage/data/d2554/l200/NJ.tre";
		String t2="/home/zhangje/storage/data/d2554/l200/TrueTree.tre";
		double dis=TreeDistance.RFdistanceScaled(t1, t2);
		System.out.println(dis);
	}

}

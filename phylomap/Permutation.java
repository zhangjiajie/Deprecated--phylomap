package phylomap;
import java.util.*;

public class Permutation {

	/**
	 * @param args
	 */
	public static int[] next(int n){
		
		Date data1=new Date();
		Random rand=new Random(data1.getTime());
		int[] a=new int[n];
		for(int i=0;i<n;i++){
			a[i]=i;
		}
		int[] b = (int[])a.clone();
		for (int k = b.length - 1; k > 0; k--) {
		    int w = rand.nextInt(k+1);
		    int temp = b[w];
		    b[w] = b[k];
		    b[k] = temp;
		}
		return b;
	}
	
	public static void printArray(int[] a) {
			for (int k = 0; k < a.length; k++)
			    System.out.print("  " + a[k]);
			System.out.println();
		    }

	
	public static void main(String[] args) {
		// TODO Auto-generated method stub
		int[] a=Permutation.next(10);
		Permutation.printArray(a);
	}

}

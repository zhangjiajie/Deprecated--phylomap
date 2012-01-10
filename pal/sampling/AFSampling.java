package pal.sampling;
import java.io.*;
import java.util.*;

import pal.alignment.*;

public class AFSampling {

	/**
	 * @param args
	 */
	
	
	public static void genAfFiles(File fin, File faln){
		try{
			ReadAlignment aln=new ReadAlignment(faln.getAbsolutePath());
			HashMap hm=new HashMap();
			for(int i=0; i<aln.getSequenceCount(); i++){
				hm.put(aln.getIdentifier(i).getName(), new Integer(i+1));
			}
			
			FileWriter fw=new FileWriter(new File(fin.getAbsolutePath()+".af"));
			FileReader fr=new FileReader(fin);
			BufferedReader bf=new BufferedReader(fr);
			String s=bf.readLine();
			Tokenizer tk;
			String t1, t2, t3;
			Integer I1, I2;
			double d3;
			double mean=0;
			double cnt=0;
			while(s!=null){
				s=s.trim();
				tk=new Tokenizer(s);
				t1=tk.nextToken();
				I1=(Integer)hm.get(t1);
				t2=tk.nextToken();
				I2=(Integer)hm.get(t2);
				t3=tk.nextToken();
				d3=Double.parseDouble(t3);
				d3=-d3;
				mean+=d3;
				cnt+=1;
				fw.write(I1+"	"+I2+"	"+d3+"\n");
				s=bf.readLine();
			}
			bf.close();
			fr.close();
			fw.close();
			System.out.println("mean of similarity: "+(mean/cnt) );
		}catch(Exception e){
			e.printStackTrace();
		}
	}
	
	public static void main(String[] args) {
		// TODO Auto-generated method stub
		File fin=new File("/home/zhangje/storage/data/NP/RAxML_distances.dis");
		File faln=new File("/home/zhangje/storage/data/NP/NP_GTR_noindel_TRUE.phy");
		//AFSampling.genAfFiles(fin, faln);
		System.out.print((Math.sqrt(2))/4);
	}

}

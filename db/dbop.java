/**
 * 
 */
package db;

import java.io.File;
import java.io.FileWriter;
import java.util.List;

/**
 * @author zhangjiajie
 *
 */
public class dbop {

	/**
	 * @param args
	 */
	
	public void toFasta(List dblist, File fout){
		
		try{
		FileWriter fw=new FileWriter(fout);
		for(int i=0;i<dblist.size();i++){
			bean b=(bean)dblist.get(i);
			if(b.getSequence().trim().length()!=0){
				String s=">"+b.getId()+"\n";
				String p=b.getSequence().trim();
				for(int k=0;k<p.length();k++){
					s+=p.substring(k, k+1);
					if(k!=0 &&(k%60)==0){
						s+="\n";
					}
				}
				s+="\n";
				s+="\n";
				fw.write(s);
			}
		}
		fw.close();
		}catch (Exception e){
			e.printStackTrace();
		}
	}
	
	public static void main(String[] args) {
		// TODO Auto-generated method stub

	}

}

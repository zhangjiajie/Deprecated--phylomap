/**
 * 
 */
package db;
import java.util.*;
import java.io.*;
import pal.alignment.*;

/**
 * @author zhangjiajie
 *
 */
public class database {

	/**
	 * @param args
	 */
	LinkedList dblist=new LinkedList();
	LinkedList tree=new LinkedList();
	ReadAlignment aln;
	int[] codebook;
	int[] clusteridx;
	
	
	public database(LinkedList dbl){
		this.dblist=dbl;
	}
	
	public database(String fphy){
		try{
			this.aln=new ReadAlignment(fphy);
			int numSeq = aln.getSequenceCount();
			for(int i=0; i<numSeq; i++){
				bean seed=new bean();
				seed.setId(aln.getIdentifier(i).getName());
				seed.setAnnotation("");
				seed.setSequence(aln.getAlignedSequenceString(i));
				dblist.add(seed);
			}
		}catch(Exception e){
			e.printStackTrace();
		}
	}
	
	public database(File fdb){
		try{
			FileReader fr=new FileReader(fdb);
			BufferedReader bf=new BufferedReader(fr);
			String s=bf.readLine();
			while(s!=null){
				String[] ss=s.split("	");
				bean seed=new bean();
				seed.setId(ss[0]);
				seed.setAnnotation(ss[1]);
				seed.setSequence(ss[2]);
				dblist.add(seed);
				s=bf.readLine();
			}
			bf.close();
			fr.close();
		}catch(Exception e){
			e.printStackTrace();
		}
	}
	
	public ReadAlignment getAlignment(){
		return this.aln;
	}
	
	public List getDB(){
		return this.dblist;
	}
	
	public void add(bean seed){
		dblist.add(seed);
	}
	
	public bean findByID(String id){
		bean br=new bean();
		for(int i=0;i<this.dblist.size();i++){
			bean b=(bean)dblist.get(i);
			if(b.getId().equals(id)){
				br=b;
				break;
			}
		}
		return br;
	}
	
	public void toFile(File fout){
		StringBuffer sb=new StringBuffer("");
		for(int i=0;i<this.dblist.size();i++){
			bean b=(bean)dblist.get(i);
			sb.append(b.toString());
		}
		try{
			FileWriter fw=new FileWriter(fout);
			fw.write(sb.toString());
			fw.close();
		}catch(Exception e){
			e.printStackTrace();
		}
	}
	
	public void assignCoord(double[][] coords){
		for(int i=0;i<dblist.size();i++){
			bean b=(bean)dblist.get(i);
			b.x=coords[0][i];
			b.y=coords[1][i];
		}
	}
	
	public void to2Dplot(File fout){
		StringBuffer sb=new StringBuffer("");
		for(int i=0;i<this.dblist.size();i++){
			bean b=(bean)dblist.get(i);
			sb.append(b.getCoords()+"\n");
		}
		
		try{
			FileWriter fw3=new FileWriter(fout);
			fw3.write(sb.toString());
			fw3.close();
		}catch(Exception e){
			e.printStackTrace();
		}
	}
	
	public void to2Dplot(File fout, String a1, String b2){
		StringBuffer sb=new StringBuffer("");
		for(int i=0;i<this.dblist.size();i++){
			bean b=(bean)dblist.get(i);
			String tag;
			if(a1.equals(b.getAnnotation())){
				tag="1";
			}else{
				tag="2";
			}
			sb.append(b.x+"	"+b.y+"	"+tag+"\n");
		}
		
		try{
			FileWriter fw3=new FileWriter(fout);
			fw3.write(sb.toString());
			fw3.close();
		}catch(Exception e){
			e.printStackTrace();
		}
	}
	
	public void to2Dplot(File fout, String a1, String b2, String c3){
		StringBuffer sb=new StringBuffer("");
		for(int i=0;i<this.dblist.size();i++){
			bean b=(bean)dblist.get(i);
			String tag;
			if(a1.equals(b.getAnnotation())){
				tag="1";
			}else if(b2.equals(b.getAnnotation())){
				tag="2";
			}else{
				tag="3";
			}
			sb.append(b.x+"	"+b.y+"	"+tag+"\n");
		}
		
		try{
			FileWriter fw3=new FileWriter(fout);
			fw3.write(sb.toString());
			fw3.close();
		}catch(Exception e){
			e.printStackTrace();
		}
	}
	
	public void setCodebook(int[] codebook){
		this.codebook=codebook;
	}
	
	public int[] getCodebook(){

		return this.codebook;
	}
	
	public void setClusteridx(int[] clusteridx){
		this.clusteridx=clusteridx;
	}
	
	public int[] getClusteridx(){
		return this.clusteridx;
	}

	public void genSampleTreeFile( File fout){
		
		for(int i=0;i<this.codebook.length;i++){
			bean b=(bean)dblist.get(codebook[i]);
			tree.add(b);
		}
		String res=this.toPhylipFormat(tree);
		try{
			FileWriter fw3=new FileWriter(fout);
			fw3.write(res);
			fw3.close();
		}catch(Exception e){
			e.printStackTrace();
		}
	}
	
	public String toPhylipFormat(List l){
		StringBuffer sb=new StringBuffer("");
		int numSpices=l.size();
		int numSites=((bean)l.get(0)).getSequence().length();
		sb.append("	"+numSpices+"	"+numSites+"\n");
		for(int i=0;i<l.size();i++){
			bean b=(bean)l.get(i);
			String acc=this.makeTo10end(b.getId());
			sb.append(acc+b.getSequence()+"\n");
		}
		
		return sb.toString();
	}
	
	public static String makeTo10end(String s){
		if(s.length()<10){
			int gap=10-s.length();
			for(int i=0;i<gap;i++){
				s+=" ";
			}
		}
		return s;
	}

	//write the leaf coords of the sampling tree
	public void genLeafCoord(File fout){
		StringBuffer sb=new StringBuffer("");
		for(int i=0;i<tree.size();i++){
			bean b=(bean)tree.get(i);
			sb.append(b.getId()+"	"+b.x+"	"+b.y+"\n");
		}
		try{
			FileWriter fw3=new FileWriter(fout);
			fw3.write(sb.toString());
			fw3.close();
		}catch(Exception e){
			e.printStackTrace();
		}
	}
	

	public static void main(String[] args) {
		// TODO Auto-generated method stub

	}

}

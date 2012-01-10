package pal.sampling;
import java.util.*;
import java.io.*;

import pal.alignment.ReadAlignment;
//import PCA.*;
//import TreeUtil.*;



public class DistacneMatrix {

	/**
	 * @param args
	 */
	private int row;
	private int col;
	private double[][] dataMatrix;
	
	public static double[][] readDistanceMatrix(File fin){
		
		try{
			FileReader fr=new FileReader(fin);
			BufferedReader bf=new BufferedReader(fr);
			String s=bf.readLine();
			int size=Integer.parseInt(s.trim());
			double[][] dm = new double[size][size];
			int cntrow=0;
			int numelement=0;
			s=bf.readLine();
			while(s!=null){
				
				int cntcol=0;
				numelement=s.length()/9;
				System.out.println(s);
				Tokenizer tker=new Tokenizer(s);
				String tk=tker.nextToken();
				tk=tker.nextToken();
				while(tk!=null){
					double ddis=Double.parseDouble(tk);
					//System.out.println(ddis);
					dm[cntrow][cntcol]=ddis;
					cntcol++;
					tk=tker.nextToken();
				}
				/*
				for(int i=1;i<numelement;i++){
					String sdis=s.substring(i*9+1,i*9+10);
					double ddis=Double.parseDouble(sdis);
					System.out.println(ddis);
					dm[cntrow][cntcol]=ddis;
					cntcol++;
				}
				*/
				
				while(cntcol<size){
					s=bf.readLine();
					tker=new Tokenizer(s);
					tk=tker.nextToken();
					while(tk!=null){
						double ddis=Double.parseDouble(tk);
						//System.out.println(ddis);
						dm[cntrow][cntcol]=ddis;
						cntcol++;
						tk=tker.nextToken();
					}
					
					/*
					numelement=s.length()/9;
					for(int i=0;i<numelement;i++){
						String sdis=s.substring(i*9,i*9+9);
						double ddis=Double.parseDouble(sdis);
						dm[cntrow][cntcol]=ddis;
						cntcol++;
					}
					*/
					
				}
				cntrow++;
				s=bf.readLine();
			}
			bf.close();
			fr.close();
			
			return dm;
		}catch(Exception e){
			e.printStackTrace();
			return null;
		}
		
		
		
		
	}

	public static double[][] readDistacneMatrixRaxml(File fdis, String falign) throws Exception{
		
		//t187 t2694 	 1.783234
		ReadAlignment aln=new ReadAlignment(falign);
		HashMap hm=new HashMap();
		for(int i=0; i<aln.getSequenceCount(); i++){
			hm.put(aln.getIdentifier(i).getName(), new Integer(i));
		}
		
		double[][] dism = new double[aln.getSequenceCount()][aln.getSequenceCount()];
		FileReader fr=new FileReader(fdis);
		BufferedReader bf=new BufferedReader(fr);
		String s=bf.readLine();
		Tokenizer tker=null;
		while(s!=null){
			tker=new Tokenizer(s);
			String sid1=tker.nextToken();
			String sid2=tker.nextToken();
			String sdis=tker.nextToken();
			//System.out.println(sid1+":"+sid2+":"+sdis);
			Integer I1=(Integer)hm.get(sid1);
			Integer I2=(Integer)hm.get(sid2);
			
			int id1=I1.intValue();
			int id2=I2.intValue();
			double dis=Double.parseDouble(sdis);
			dism[id1][id2]=dis;
			dism[id2][id1]=dis;
			s=bf.readLine();
		}
		
		bf.close();
		fr.close();
		return dism;
	}
	
	public static double[][] readDistacneMatrixRaxml(File fdis, ReadAlignment aln) throws Exception{
		
		//t187 t2694 	 1.783234
		//ReadAlignment aln=new ReadAlignment(falign);
		HashMap hm=new HashMap();
		for(int i=0; i<aln.getSequenceCount(); i++){
			hm.put(aln.getIdentifier(i).getName(), new Integer(i));
		}
		
		double[][] dism = new double[aln.getSequenceCount()][aln.getSequenceCount()];
		FileReader fr=new FileReader(fdis);
		BufferedReader bf=new BufferedReader(fr);
		String s=bf.readLine();
		Tokenizer tker=null;
		while(s!=null){
			tker=new Tokenizer(s);
			String sid1=tker.nextToken();
			String sid2=tker.nextToken();
			String sdis=tker.nextToken();
			//System.out.println(sid1+":"+sid2+":"+sdis);
			Integer I1=(Integer)hm.get(sid1);
			Integer I2=(Integer)hm.get(sid2);
			
			int id1=I1.intValue();
			int id2=I2.intValue();
			double dis=Double.parseDouble(sdis);
			dism[id1][id2]=dis;
			dism[id2][id1]=dis;
			s=bf.readLine();
		}
		
		bf.close();
		fr.close();
		return dism;
	}
	
	public void compare(double[][] dm1, double[][] dm2, File fout){
		StringBuffer sb=new StringBuffer("");
		int size=dm1.length;
		int size2=dm2.length;
		if(size!=size2){
			System.out.println("Matrix not the same dimension!!!!!!!!!!!!!!!!!!!!!!!!!!");
		}else{
			for(int i=0;i<size;i++){
				for(int j=i;j<size;j++){
					sb.append(dm1[i][j]+"	"+dm2[i][j]+"\n");
				}
			}
		}
		try{
			FileWriter fw3=new FileWriter(fout,true);
			fw3.write(sb.toString());
			fw3.close();
		}catch(Exception e){
			e.printStackTrace();
		}
		
	}
	
	public void toFile(String fout){
		this.row=dataMatrix.length;
		this.col=dataMatrix[0].length;
		try{
			FileWriter fw3=new FileWriter(new File (fout),true);
			fw3.write(this.row+"	"+this.col+"\n");
			for(int i=0;i<row;i++){
				StringBuffer sb=new StringBuffer("");
				for (int j=0;j<col;j++){
					sb.append(dataMatrix[i][j]+"	");
				}
				fw3.write(sb.toString().trim()+"\n");
			}
			fw3.close();
		}catch(Exception e){
			e.printStackTrace();
		}
	}
	
	public void readFromFile (String fin){
		try{
			FileReader fr=new FileReader(new File(fin));
			BufferedReader bf=new BufferedReader(fr);
			String s=bf.readLine();
			String[] ss=s.split("\\s");
			this.row=Integer.parseInt(ss[0]);
			this.col=Integer.parseInt(ss[1]);
			this.dataMatrix=new double[this.row][this.col];
			int rowcnt=0;
			s=bf.readLine();
			while(s!=null){
				s=s.trim();
				ss=s.split("\\s");
				for(int i=0; i<ss.length; i++){
					dataMatrix[rowcnt][i]=Double.parseDouble(ss[i]);
				}
				rowcnt++;
				s=bf.readLine();
			}
			bf.close();
			fr.close();
		}catch(Exception e){
			e.printStackTrace();
		}
	}
	
	public List readPhylip_seqeunce_alignemnt(File fin){
		LinkedList seqs=new LinkedList();
		
		try{
			FileReader fr=new FileReader(fin);
			BufferedReader bf=new BufferedReader(fr);
			String s=bf.readLine();
			s=bf.readLine();
			int count=0;
			String seq="";
			String acc="";
			while(s!=null){
				seq=s.substring(10);
				seqs.add(seq);
				s=bf.readLine();
			}
			
			bf.close();
			fr.close();
			
			
		}catch(Exception e){
			e.printStackTrace();
		}
		
		return seqs;
		
	}
	
	public void delete_Gap_phylip(File fin, File fout){
		StringBuffer sb=new StringBuffer();
		try{
			FileReader fr=new FileReader(fin);
			BufferedReader bf=new BufferedReader(fr);
			String s=bf.readLine();
			sb.append(s+"\n");
			s=bf.readLine();
			
			int count=0;
			String seq="";
			String acc="";
			while(s!=null){
				acc=s.substring(0, 10);
				seq=s.substring(35,168)+s.substring(169, 771);
				//seqs.add(seq);
				sb.append(acc+seq+"\n");
				s=bf.readLine();
			}
			
			bf.close();
			fr.close();
			
			FileWriter fw3=new FileWriter(fout);
			fw3.write(sb.toString());
			fw3.close();
		}catch(Exception e){
			e.printStackTrace();
		}
	}
	
	/*
	public void compareMetric(){
		PCoA pcoa=new PCoA();
		File fin_matrix_phy=new File("/storage/disk3/test/outfile");
		File fin_seqs_phy=new File("/storage/disk3/test/sample100_nogap.phy");
		File fout=new File("/storage/disk3/test/metric_compare3.dat");
		double[][] dm_phy=this.readDistanceMatrix(fin_matrix_phy);
		double[][] dm_pdis=pcoa.dataMatrix_p(this.readPhylip_seqeunce_alignemnt(fin_seqs_phy));
		this.compare(dm_pdis, dm_phy, fout);
		
	}
	*/
	
	//scan one sequence for possible kwords
	private HashSet<String> scanOneSeq(HashSet<String> hs, String seq, int len){
		seq=new String (seq.replaceAll("-", ""));
		for(int i=0;i<seq.length()-len+1;i++){
			hs.add(seq.substring(i, i+len));
		}
		return hs;
	}
	
	private HashSet<String> scanSeqSet(HashSet<String> hs, List<String> seqlist, int len){
		for(int i=0; i<seqlist.size(); i++){
			hs=this.scanOneSeq(hs, seqlist.get(i), len);
		}
		return hs;
	}
	
	private double[] seq2vec (String seq, String[] kmer, int k){
		seq=new String (seq.replaceAll("-", ""));
		int size = kmer.length;
		int len = seq.length();
		int denom = len-k+1;
		double[] vec = new double[size];
		String[] skmer=new String[denom];
		for(int i=0; i<denom; i++){
			skmer[i]=seq.substring(i, i+k);
		}
		String currkmer;
		for(int i=0; i<size; i++){
			double cnt=0;
			currkmer=kmer[i];
			for(int j=0;j<denom;j++){
				if(currkmer.equals(skmer[j])){
					cnt+=1.0;
				}
			}
			vec[i]=cnt/((double)(denom));
		}
		return vec;
	}

	public double[][] genCVMatrix(ReadAlignment aln, int k){
		HashSet hs =new HashSet();
		String[] kmer;
		int numSeq=aln.getSequenceCount();
		int numSite=aln.getSiteCount();
		for(int i=0; i<numSeq; i++){
			String seq=aln.getAlignedSequenceString(i);
			hs=this.scanOneSeq(hs, seq, k);
		}
		Object[] objs= hs.toArray();
		kmer=new String[objs.length];
		for(int i=0; i<hs.size(); i++){
			kmer[i]=(String)objs[i];
		}
	
		//this.kmer=(String[])this.hs.toArray();
		int size=kmer.length;
		this.dataMatrix=new double[numSeq][size];
		for(int i=0; i<numSeq; i++){
			String seq=aln.getAlignedSequenceString(i);
			this.dataMatrix[i]=this.seq2vec(seq, kmer, k);
		}
		
		return this.dataMatrix;
	}
	
	public static void main(String[] args) {
		// TODO Auto-generated method stub
		//File fin=new File("/storage/disk2/Output/outfile");
		DistacneMatrix dm=new DistacneMatrix();
		//dm.compareMetric();
		
		//File fin_seqs_phy=new File("/storage/disk3/test/sample100.phy");
		//File fout=new File("/storage/disk3/test/sample100_nogap.phy");
		//dm.delete_Gap_phylip(fin_seqs_phy, fout);
		//File fout=new File("/storage/disk2/Output/padiscomparesmith");
		//dm.toFile(dm.readDistanceMatrix(fin), fout);
		//fastaDB fd=new fastaDB();
		//dm.compare(dm.readDistanceMatrix(fin), fd.getDisMatrix("PA"), fout);
		
		String s="t187 t2694 	 1.783234";
		String[] ss=s.split("\\s");
		for(int i=0; i<ss.length; i++){
			System.out.println(ss[i]);
		}
	
	}

}

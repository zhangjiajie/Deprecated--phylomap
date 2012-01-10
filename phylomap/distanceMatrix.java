package phylomap;

import java.util.*;
import java.io.*;
import java.sql.*;

import pal.alignment.*;
import pal.sampling.Tokenizer;
import pal.substmodel.*;
import pal.datatype.*;


public class distanceMatrix {

	/**
	 * @param args
	 */
	private int row;
	private int col;
	private double[][] dataMatrix;
	
	public double[][] readDistanceMatrix(File fin){
		
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
				s=s.trim();
				if(!"".equals(s)){
					int cntcol=0;
					numelement=s.length()/10;
					//System.out.println(numelement);
					for(int i=1;i<numelement;i++){
						String sdis=s.substring(i*10,i*10+10);
						double ddis=Double.parseDouble(sdis);
						//System.out.println(ddis);
						dm[cntrow][cntcol]=ddis;
						cntcol++;
					}
					
					while(cntcol<size){
						s=bf.readLine();
						
						numelement=s.length()/10;
						for(int i=0;i<numelement;i++){
							String sdis=s.substring(i*10,i*10+10);
							double ddis=Double.parseDouble(sdis);
							dm[cntrow][cntcol]=ddis;
							cntcol++;
						}
						//System.out.println("should have no id: "+s+"	"+cntcol);
					}
					cntrow++;
				}
				s=bf.readLine();
			}
			bf.close();
			fr.close();
			
			this.dataMatrix = dm;
			
			return dm;
		}catch(Exception e){
			e.printStackTrace();
			return null;
		}
		
	}

	public double[][] readPhylipDistanceMatrix(File fin){
		
		try{
			FileReader fr=new FileReader(fin);
			BufferedReader bf=new BufferedReader(fr);
			String s=bf.readLine();
			int size=Integer.parseInt(s.trim());
			double[][] dm = new double[size][size];
			int cntrow=0;
			s=bf.readLine();
			while(s!=null){
				
				int cntcol=0;
				System.out.println(s);
				Tokenizer tker=new Tokenizer(s);
				String tk=tker.nextToken();
				tk=tker.nextToken();
				while(tk!=null){
					double ddis=Double.parseDouble(tk);
					dm[cntrow][cntcol]=ddis;
					cntcol++;
					tk=tker.nextToken();
				}
				
				while(cntcol<size){
					s=bf.readLine();
					tker=new Tokenizer(s);
					tk=tker.nextToken();
					while(tk!=null){
						double ddis=Double.parseDouble(tk);
						dm[cntrow][cntcol]=ddis;
						cntcol++;
						tk=tker.nextToken();
					}
				}
				cntrow++;
				s=bf.readLine();
			}
			bf.close();
			fr.close();
			this.dataMatrix = dm;
			return dm;
		}catch(Exception e){
			e.printStackTrace();
			return null;
		}
	}
	
	public double[][] readRAxMLDistacneMatrix(File fdis, ReadAlignment aln) throws Exception{
		
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
		this.dataMatrix = dism;
		return dism;
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
					this.dataMatrix[rowcnt][i]=Double.parseDouble(ss[i]);
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
	
	public void compare(double[][] dm1, double[][] dm2, File fout){
		StringBuffer sb=new StringBuffer("");
		int size=dm1.length;
		int size2=dm2.length;
		if(size!=size2){
			System.out.println("Matrix not the same dimension!");
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
	
	/*
	public void toFile(double[][] dm, File fout){
		
		try{
			FileWriter fw3=new FileWriter(fout,true);
			int k=dm.length;
			for(int i=0;i<k;i++){
				StringBuffer sb=new StringBuffer("");
				for (int j=0;j<k;j++){
					sb.append(dm[i][j]+"	");
				}
				fw3.write(sb.toString().trim()+"\n");
			}
			fw3.close();
		}catch(Exception e){
			e.printStackTrace();
		}
	}
	*/
	
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
	
	/*
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
			seq = seq.replaceAll("-", "");
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
			seq = seq.replaceAll("-", "");
			this.dataMatrix[i]=this.seq2vec(seq, kmer, k);
		}
		
		return this.dataMatrix;
	}

	
	double[] A= {0.3535533905932738, 0.5, 0.0};
	double[] G= {0.3535533905932738, -0.5, 0.0};
	double[] T=	{-0.3535533905932738, 0.0, 0.5};
	double[] U=	{-0.3535533905932738, 0.0, 0.5};
	double[] C=	{-0.3535533905932738, 0.0, -0.5};
	double[] Space= {0.0, 0.0, 0.0};
	HashMap<String, double[]> hm = new HashMap<String, double[]>();
	
	public double[] translate(String seq){
		int len=seq.trim().length()*3;
		double[] vec=new double[len];
		int size = seq.trim().length();
		String nn;
		int cnt=0;
		double[] v;
		for(int i=0; i<size; i++){
			nn=seq.substring(i,i+1);
			v= hm.get(nn);
			if(v==null){
				v=Space;
			}
			for(int k=0;k<v.length;k++){
				vec[cnt]=v[k];
				cnt++;
			}
		}
		return vec;
	}
	
	public double[][] genSVMatrix(ReadAlignment aln){
		hm.put("A", A);
		hm.put("a", A);
		hm.put("G", G);
		hm.put("g", G);
		hm.put("T", T);
		hm.put("t", T);
		hm.put("U", U);
		hm.put("u", U);
		hm.put("C", C);
		hm.put("c", C);
		hm.put("-", Space);
		int numSeq=aln.getSequenceCount();
		int numSite=aln.getSiteCount();
		this.dataMatrix = new double[numSeq][numSite*3];
		for(int i=0; i<numSeq; i++){
			String seq=aln.getAlignedSequenceString(i);
			double[] vec=this.translate(seq);
			this.dataMatrix[i]=vec;
			//System.out.println("Translating sequence: "+i);
		}
		return this.dataMatrix;

}
	
	public double norm2vector(double[] v1, double[] v2){
		double dis=0;
		for(int i=0;i<v1.length;i++){
			dis=dis+(v1[i]-v2[i])*(v1[i]-v2[i]);
		}
		return Math.sqrt(dis);
	}
	
	public double[][] genNbyNDistanceMatrix(double[][] dataMatrix){
		
		int n = dataMatrix.length;
		double[][] dism = new double[n][n];
		for(int i=0; i<n; i++){
			for(int j=i; j<n; j++){
				if(i == j){
					dism[i][j] = 0;
				}else{
					dism[i][j] = this.norm2vector(dataMatrix[i], dataMatrix[j]);
					dism[j][i] = dism[i][j];
				}
			}
		}
		
		return dism;
	}
	
	public static void main(String[] args) {
		// TODO Auto-generated method stub
		//PushbackReader in;
		//in = FileIO.openIn(indataName);
		String fin = "/Users/zhangjiajie/Desktop/sample.phy";
		Alignment raw = null;
		try
		{
			raw = new ReadAlignment(fin);
		}
		catch (AlignmentParseException e)
		{
			System.out.println("Error: Alignment parsing problem");
			System.exit(1);
		}
		catch (IOException e)
		{
			System.out.println("Error: File not found (IO error)");
			System.exit(1);
		}	
		//FileIO.close(in);
		System.out.println("Contains " + raw.getSequenceCount() + " sequences of length " + raw.getSiteCount());
		System.out.println("Likely content: " + raw.getDataType().getDescription() + " data");
		
		Options options = new Options();
		options.setMLDIST();
		options.dtyp = raw.getDataType().getTypeID();
		
		if(options.dtyp == 1){
			options.smodel = AminoAcidModelID.WAG;
		}else{
			options.smodel = NucleotideModelID.GTR;
		}
		
		// Start computation
		java.util.Date date=new java.util.Date();
		Timestamp timeStamp = new Timestamp(date.getTime());

		// Compute distances
		//raw.setDataType( DataType.Utils.getInstance(options.dtyp) );
		
		SitePattern sitePattern = new SitePattern(raw);
		
		pal.distance.DistanceMatrix mat = null;
		
		System.out.println();
		System.out.println("Computing distance matrix");
		
		SubstitutionModel model = null;
	
		RateMatrix rmat = RateMatrixUtils.getInstance(options.dtyp, options.smodel, options.params, AlignmentUtils.estimateFrequencies(raw));
	
		RateDistribution rdist = null;
		if (options.rmodel == 0) rdist = new UniformRate();
		if (options.rmodel == 1) rdist = new GammaRates(options.alphaCats, options.alpha);	
		if (options.rmodel == 2) rdist = new InvariableSites(options.fracInv);	
		//model = new SubstitutionModel(rmat, rdist);	

		
	}

}


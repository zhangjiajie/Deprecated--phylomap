import java.util.*;
import java.io.*;

import db.*;
import phylomap.*;
import mds.*;

public class phylomaps1 {

	/**
	 * @param args
	 */
	private static int atoi(String s)
	{
		return Integer.parseInt(s);
	}

	public void  analysisA(File fdb, File fMatrix, File f2dplot, File eigenvalueout, File fsample, File fsample_coord, int numsample){
		
		/*Input
		 * fdb: the database, with each line a recored, each record should have 3 fields separated by TAB, first field is ID with less than 10 char, 
		 * second is annotation; third is the sequence 
		 * fMatrix: calculated by the phylip package program: protdist
		 * numsample: number of sequence centers
		 *  
		 *Output:  
		 *f2dplot: file for plot in matlab  
		 *eigenvalueout: eigenvalues from PCoA  
		 *fsample: the sequence centers find by clustering in fasta format, this file will be used by phylip to build a sample tree (NJ tree only at the moment)
		 *fsample_coord: file for plot in matlab
		 *
		 */
		
		System.out.println("Input database:"+fdb.getAbsolutePath() );
		System.out.println("Distance matrix:"+fMatrix.getAbsolutePath());
		
		//init db
		database db=new database(fdb);
		
		//read distance matrix 
		distanceMatrix dm=new distanceMatrix();
		double[][] dataMatrix=dm.readDistanceMatrix(fMatrix);
		
		System.out.println("Start PCoA ..........");
		//PCoA
		PCoA pca=new PCoA();
		//get first 2 axis
		double[][] coords=pca.PcoordPhylip(dataMatrix);
		//assign each bean a coordinate
		db.assignCoord(coords);
		//save to file
		db.to2Dplot(f2dplot);
		//save eigenvalues
		double eigen1=pca.accumulateEigenvaluesTofile(eigenvalueout);
		
		System.out.println("Start NeuralGas Clustering ..........");
		//Vector quantization
		NeuralGas ng=new NeuralGas(dataMatrix, numsample, 100);
		ng.learn();
		int[] codebook=ng.getCodebookvector();
		db.setCodebook(codebook);
		db.setClusteridx(ng.getClusterindex());
		db.genSampleTreeFile(fsample);
		db.genLeafCoord(fsample_coord);
	
	}
	
	private static void exit_with_help()
	{
		System.out.print(
		 "Usage: phylomaps1 -d database -m matrix -n num_cluster_center\n"
		+"-d database file name\n"
		+"-m distance matrix file name\n"
		+"-n number of sampling sequences (default 30)\n"
		);
		System.exit(1);
	}
	
	
	private void run(String argv[]) throws IOException
	{		
		String sfdb="";
		String fmatrix="";
		String matrixType="raxml";
		int numCenter=30;
		// parse options
		for( int i=0;i<argv.length;i++)
		{
			if(argv[i].charAt(0) != '-') break;
			if(++i>=argv.length)
				exit_with_help();
			switch(argv[i-1].charAt(1))
			{
				
				case 'd':
					sfdb = argv[i];
					break;
				case 'm':
					fmatrix= argv[i];
					break;
				case 'P':
					matrixType = "phylip";
					break;
				case 'X':
					matrixType = "raxml";
					break;
				case 'n':
					numCenter = atoi(argv[i]);
					break;
				default:
					System.err.print("Unknown option: " + argv[i-1] + "\n");
					exit_with_help();
			}
		}
		
		if("".equals(sfdb) || "".equals(fmatrix)){
			exit_with_help();
		}else{
		
			String sf2dplot=sfdb+".2dplot";
			String seigenvalueout=sfdb+".eigen";
			String sfsample=sfdb+".sample";
			String sfnamecoord=sfdb+".namecoord";
			String sfsample_coord=sfdb+".samplecoord";
			
			//input
			File fdb=new File(sfdb);
			
			File fMatrix=new File(fmatrix);
			
			//output
			File fsample_coord=new File(sfsample_coord);
			File f2dplot=new File(sf2dplot);
			File eigenvalueout=new File(seigenvalueout);
			File fsample=new File(sfsample);
			File fnamecoord=new File(sfnamecoord);
			
			//run
			this.analysisA(fdb, fMatrix, f2dplot, eigenvalueout, fsample, fsample_coord, numCenter);
		}
	
	}
	public static void main(String[] args) throws IOException {
		// TODO Auto-generated method stub
		phylomaps1 pms1=new phylomaps1();
		pms1.run(args);
	}

}

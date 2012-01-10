import java.util.*;
import java.io.*;

import db.*;
import pal.alignment.SubAlignment;
import pal.distance.DistanceMatrix;
import pal.misc.Identifier;
import pal.misc.SimpleIdGroup;
import pal.sampling.Error;
import pal.tree.ExtractSubtree;
import pal.tree.NeighborJoiningTree;
import phylomap.*;
import mds.*;

public class phylomap {

	/**
	 * @param args
	 */
	private static int atoi(String s)
	{
		return Integer.parseInt(s);
	}

	public void analysis(database db, double[][] dataMatrix, String outputname, File f2dplot, File eigenvalueout, File fsample, File fsample_coord, int numsample, String clustering_method, boolean isDistanceMatrix){
		
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
		
		System.out.println("Start PCoA ..........");
		//PCoA
		PCoA pca=new PCoA();
		//get first 2 axis
		double[][] coords = null;
		double[][] distanceMatrix = null;
		if(isDistanceMatrix){
			coords=pca.PcoordPhylip(dataMatrix);
		}else{
			distanceMatrix dm = new distanceMatrix();
			System.out.println("Calculate distance matrix for PCoA ..........");
			distanceMatrix = dm.genNbyNDistanceMatrix(dataMatrix);
			coords = pca.PcoordPhylip(distanceMatrix);
		}
		
		//assign each bean a coordinate
		db.assignCoord(coords);
		//save to file
		db.to2Dplot(f2dplot);
		//save eigenvalues
		double eigen1=pca.accumulateEigenvaluesTofile(eigenvalueout);
		System.out.println("PCoA end");
		
		System.out.println("Start "+clustering_method+" Clustering ..........");
		int[] codebook = null;
		if("NeuralGas".equals(clustering_method)){
			NeuralGas ng=new NeuralGas(dataMatrix, numsample, 100);
			ng.learn();
			codebook=ng.getCodebookvector();
			db.setCodebook(codebook);
			db.setClusteridx(ng.getClusterindex());
			db.genSampleTreeFile(fsample);
			db.genLeafCoord(fsample_coord);
		}else if("PAM".equals(clustering_method)){
			PAM pam = new PAM(dataMatrix, numsample);
			codebook = pam.getCodebook();
			db.setCodebook(codebook);
			db.setClusteridx(pam.getAssignment());
			db.genSampleTreeFile(fsample);
			db.genLeafCoord(fsample_coord);
		}else{
			System.err.print("Unkown clustering method!" + "\n");
			exit_with_help();
		}
		
		
		//Vector quantization
		
		int numSeqs=numsample;
		int numSites=db.getAlignment().getSiteCount();
		
		// Reserve memory
		SimpleIdGroup idGroup = new SimpleIdGroup(numSeqs);
		char[][] data = new char[numSeqs][numSites];
		//HashSet hs=new HashSet();
		
		System.out.println("Selected sequences index:");
		for(int i=0; i<numsample; i++){
			System.out.print(codebook[i]+",");
			data[i]=db.getAlignment().getAlignedSequenceString(codebook[i]).toCharArray();
			Identifier id=db.getAlignment().getIdentifier(codebook[i]);
			//hs.add(id.getName());
			idGroup.setIdentifier(i, id);
		}
		System.out.println();
		SubAlignment subaln=new SubAlignment(numSeqs, numSites, data, idGroup);
		//Error err=new Error(dataMatrix, this.sampleAlignmentFile.getAbsolutePath(), db.getAlignment());
		double[][] sampleDataMatrix=new double[subaln.getSequenceCount()][subaln.getSequenceCount()];
		
		int[] idx = new int[subaln.getSequenceCount()];
		for(int i=0; i<subaln.getSequenceCount(); i++){
			idx[i]=db.getAlignment().whichIdNumber(subaln.getIdentifier(i).getName());
		}
		
		if(isDistanceMatrix){
			for(int i=0; i<subaln.getSequenceCount(); i++){
				for (int j=0; j<subaln.getSequenceCount(); j++){
					sampleDataMatrix[i][j]=dataMatrix[idx[i]][idx[j]];
				}
			}
		}else{
			for(int i=0; i<subaln.getSequenceCount(); i++){
				for (int j=0; j<subaln.getSequenceCount(); j++){
					sampleDataMatrix[i][j]=distanceMatrix[idx[i]][idx[j]];
				}
			}
		}

		
		System.out.println(clustering_method+" Clustering end");
		
		System.out.println("Building NJ tree ......");
		
		DistanceMatrix dism = new DistanceMatrix(sampleDataMatrix, subaln.getIdGroup());
	
		NeighborJoiningTree njt = new NeighborJoiningTree (dism);
		ExtractSubtree.writeSubTree(new File(outputname+".tre"), njt);
		
		System.out.println("Start mapping the sampling tree onto the PCoA results ..........");
		pnp Pnp=new pnp();
		try{
			HashMap Name_Coords=Pnp.runPnp(db, new File(outputname+".namecoord"), new File(outputname+".tre"), new File(outputname+".samplecoord"), numsample, -1, 1, -1, 1, 10000, new File(outputname+".cc"), new File(outputname+".beziercc"), 1);
			Pnp.Name_Coords_toString(Name_Coords, new File(outputname+".innercoord"));
			postOrder po=new postOrder();
			node root=Pnp.tree_root;
			po.genTreePlot(new File(outputname+".treecoord"), root, Name_Coords, Pnp.hm_controlpoint);
		}catch(Exception e){
			e.printStackTrace();
		}
		
	}
	
	private static void exit_with_help()
	{
		System.out.print(
		 "Usage: phylomap -d alignment -m matrix -x raxml -n num_cluster_center -k kmer_length\n"
		+"-d alignment\n"
		+"-m data matrix file name (this file is needed only when using data matrix type 'phylip' or 'raxml')\n"
		+"-x data matrix type: 'phylip' or 'raxml' or 'CV' or 'SV'\n"
		+"-c clustering method: 'NeuralGas' or 'PAM', when using data matrix type: 'phylip' or 'raxml', an alternative clustering method PAM can be used for better results.\n"
		+"   Note PAM is extremly slow for large number of sampling sequences on large data set, and can only be used with data matrix type: 'phylip' or 'raxml', the default option is 'NeuralGas'\n"
		+"-n number of sampling sequences (default 30)\n"
		+"-k length of kmer string if use CV data matrix type (default 5)\n"
		
		+"Please contact bestzhangjiajie@gmail.com if you have any questions.\n"
		);
		System.exit(1);
	}
	
	private void run(String argv[]) throws IOException
	{		
		String sfdb="";
		String fmatrix="";
		String matrixType="raxml"; //raxml, phylip, CV, SV, random
		String clusterMethod="NeuralGas"; //NeuralGas, PAM
		int numCenter=30;
		int k=5;
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
				case 'x':
					matrixType = argv[i];
					break;
				case 'c':
					clusterMethod = argv[i];
					break;	
				case 'n':
					numCenter = atoi(argv[i]);
					break;
				case 'k':
					k = atoi(argv[i]);
					break;	
				default:
					System.err.print("Unknown option: " + argv[i-1] + "\n");
					exit_with_help();
			}
		}
		
		//check options
		if("".equals(sfdb)){
			exit_with_help();
		}else{
			String sf2dplot=sfdb+".2dplot";
			String seigenvalueout=sfdb+".eigen";
			String sfsample=sfdb+".sample";
			String sfnamecoord=sfdb+".namecoord";
			String sfsample_coord=sfdb+".samplecoord";
			boolean isDistanceMatrix = false;
			
			//input
			File fdb=new File(sfdb);
			
			File fMatrix=new File(fmatrix);
			
			//output
			File fsample_coord=new File(sfsample_coord);
			File f2dplot=new File(sf2dplot);
			File eigenvalueout=new File(seigenvalueout);
			File fsample=new File(sfsample);
			File fnamecoord=new File(sfnamecoord);
			
			//print out all parameters
			System.out.println("Alignment: " + sfdb);
			System.out.println("Distance Matrix: " + fmatrix);
			System.out.println("Distance Matrix type: " + matrixType);
			System.out.println("Clustering method: " + clusterMethod);
			System.out.println("Number of sampling sequences: " + numCenter);
			if ("CV".equals(matrixType)){System.out.println("Kmer length: " + k);}
			
			database db = null;
			double[][] dataMatrix = null;
			if("phylip".equals(matrixType)){
				if("".equals(fmatrix)){
					exit_with_help();
				}else{
					db = new database(sfdb);
					distanceMatrix dism = new distanceMatrix();
					dataMatrix = dism.readPhylipDistanceMatrix(fMatrix);
					isDistanceMatrix = true;
				}
			}else if("raxml".equals(matrixType)){
				if("".equals(fmatrix)){
					exit_with_help();
				}else{
					db = new database(sfdb);
					distanceMatrix dism = new distanceMatrix();
					try{
						dataMatrix = dism.readRAxMLDistacneMatrix(fMatrix, db.getAlignment());
					}catch(Exception e){
						e.printStackTrace();
					}
					isDistanceMatrix = true;
				}
			}else if("CV".equals(matrixType)){
				if("NeuralGas".equals(clusterMethod)){
					db = new database(sfdb);
					distanceMatrix dism = new distanceMatrix();
					try{
						System.out.println("Transforming the sequences into feature frequences, this can take a while if -k is larger than 7 ......" + sfdb);
						dataMatrix = dism.genCVMatrix(db.getAlignment(), k);
					}catch(Exception e){
						e.printStackTrace();
					}
					isDistanceMatrix = false;
				}else{
					System.err.print("CV data matrix can only use NeuralGas as clustering method" + "\n");
					exit_with_help();
				}
				
			}else if("SV".equals(matrixType)){
				if("NeuralGas".equals(clusterMethod)){
					db = new database(sfdb);
					distanceMatrix dism = new distanceMatrix();
					try{
						dataMatrix = dism.genSVMatrix(db.getAlignment());
					}catch(Exception e){
						e.printStackTrace();
					}
					isDistanceMatrix = false;
				}else{
					System.err.print("SV data matrix can only use NeuralGas as clustering method" + "\n");
					exit_with_help();
				}
			}else if("random".equals(matrixType)){
				System.err.print("Not supported so far!" + "\n");
				exit_with_help();
				
			}else{
				System.err.print("Unkown distance matrix type!" + "\n");
				exit_with_help();
			}
			
			
			//run
			this.analysis(db, dataMatrix, sfdb, f2dplot, eigenvalueout, fsample, fsample_coord, numCenter, clusterMethod, isDistanceMatrix);
			
			System.out.println();
			System.out.println("The sampling alignment saved to:	" +sfdb+".sample");
			System.out.println("The sampling NJ tree saved to:	" +sfdb+".tre");
			System.out.println("The accumulated eigenvalue saved to:	" +sfdb+".eigen");
			System.out.println("The input alignment first 2 PCoA axis coordinates saved to:	" +sfdb+".2dplot");
			System.out.println("The sampling first 2 PCoA axis coordinates saved to:	" +sfdb+".samplecoord");
			
			System.out.println("The mapped tree inner nodes coordinates saved to:	" +sfdb+".innercoord");
			System.out.println("The mapped tree taxa name coordinates saved to:	" +sfdb+".namecoord");			
			System.out.println("The bezier coordinates for drawing the mapped tree saved to:	" +sfdb+".treecoord");			
			System.out.println("A side by side compare of branch length after mapping saved to:	" +sfdb+".cc");
			System.out.println("A side by side compare of branch length after mapping with bezier compensation saved to:	" +sfdb+".beziercc");
			
			System.out.println();
			System.out.println("Please contact bestzhangjiajie@gmail.com if you have any questions.");
		}
		
	}
	
	public static void main(String[] args) throws IOException {
		// TODO Auto-generated method stub
		phylomap pms=new phylomap();
		pms.run(args);
	}

}

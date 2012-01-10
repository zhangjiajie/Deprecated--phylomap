import java.util.*;
import java.io.*;

import db.*;
import phylomap.*;
import mds.*;

public class phylomaps2 {

	/**
	 * @param args
	 */
	private static int atoi(String s)
	{
		return Integer.parseInt(s);
	}
	
	private static void exit_with_help()
	{
		System.out.print(
		 "Usage: phylomaps2 -d database -c *.namecoord -n num_cluster_center -t tree -k *.samplecoord\n"
		+"-d database file name\n"
		+"-c annotation coordinate, generated by phylomaps1: *.namecoord\n"
		+"-n number of sampling sequences (default 30)\n"
		+"-t unrooted NJ tree generated by Phylip (do not support other types of tree now)\n"
		+"-k sampling sequence 2d coordinates generated by phylomaps1: *.samplecoord\n"
		+"-i  max iters for tree mapping (default 10000)\n"
		);
		System.exit(1);
	}
	
	
	private void run(String argv[]){
		String sfdb="";
		String sfnamecoord="";
		String streefile="";
		String sfsample_coord="";
		int numCenter=30;
		int maxiters=10000;
		
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
				case 'c':
					sfnamecoord= argv[i];
					break;
				case 'n':
					numCenter = atoi(argv[i]);
					break;
				case 't':
					streefile = argv[i];
					break;
				case 'k':
					sfsample_coord = argv[i];
					break;
				case 'i':
					maxiters = atoi(argv[i]);
					break;	
				default:
					System.err.print("Unknown option: " + argv[i-1] + "\n");
					exit_with_help();
			}
		}
		
		if("".equals(sfdb) || "".equals(sfnamecoord) || "".equals(streefile) || "".equals(sfsample_coord)){
			exit_with_help();
		}else{
		
			File fdb=new File(sfdb);
			File fnamecoord=new File(sfnamecoord);
			File treefile=new File(streefile);
			File fsample_coord=new File(sfsample_coord);
			File foutcc=new File(sfdb+".cc");
			File foutcc_bezier=new File(sfdb+".beziercc");
			File foutinnercoord=new File(sfdb+".innercoord");
			File fouttreecoord=new File(sfdb+".treecoord");
			
			pnp Pnp=new pnp();
			try{
				HashMap Name_Coords=Pnp.runPnp(fdb, fnamecoord, treefile, fsample_coord, numCenter, -1, 1, -1, 1, 10000, foutcc, foutcc_bezier, 1);
				Pnp.Name_Coords_toString(Name_Coords, foutinnercoord);
				postOrder po=new postOrder();
				node root=Pnp.tree_root;
				po.genTreePlot(fouttreecoord, root, Name_Coords, Pnp.hm_controlpoint);
			}catch(Exception e){
				e.printStackTrace();
			}
			
		}
		
	}
	
	
	public static void main(String[] args) {
		// TODO Auto-generated method stub
		phylomaps2 pms2=new phylomaps2();
		pms2.run(args);
	}

}

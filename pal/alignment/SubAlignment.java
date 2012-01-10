package pal.alignment;
import java.util.*;
import java.io.*;
import pal.misc.*;

public class SubAlignment extends ReadAlignment {

	/**
	 * @param args
	 */
	public SubAlignment(int numSeq, int numSite, char[][] data, SimpleIdGroup idGroup){
		this.numSeqs=numSeq;
		this.numSites=numSite;
		this.data=data;
		this.idGroup = idGroup;
	}
	
	public void toFile(File fout){
		try{
			PrintWriter pw=new PrintWriter(fout);
			AlignmentUtils.printInterleaved(this, pw);
			pw.close();
		}catch(Exception e){
			e.printStackTrace();
		}
	}
	
	public static void main(String[] args) {
		// TODO Auto-generated method stub

	}

}

// Options.java
//
// (c) 1999-2001 Korbinian Strimmer
//
// This package may be distributed under the
// terms of the GNU General Public License (GPL)


package phylomap;

import pal.datatype.*;
import pal.substmodel.*;
import pal.io.*;

import java.io.*;

/**
 * used by MLTREE, MLDIST, EVOLVE, DISTTREE, REWRITE
 * to display and set options
 *
 * @version $Id:$
 *
 * @author Korbinian Strimmer
 */
public class Options
{
	//
	// Public stuff
	//

	public Options()
	{
		dtyp = DataType.NUCLEOTIDES;
		smodel = NucleotideModelID.F84; //F84
		rmodel = 0; //Uniform rate
		alphaCats = 8;
		alpha = 1.0;
		fracInv = 0.0;
		treeTest = 0;
		
		numDataSets = 1;
		numSites = 1000;
		format = 0;
		distMethod = 0;
		treeMethod = 0;
		numBootstraps = 0;
		branchConstraint = 0;
		dropSites = 0;
		
		params = new double[5];
		resetParams();
		
		freq = null;
		userFreqs = false;
		optTree = true;
		optModel = false;
		useModel = true;
		siteLik = false;
		weightedLS = true;
		dropGaps = true;
		jumble = false;
		bootstrap = false;
	}

	//
	// Friendly stuff
	//

	// Parameters
	int dtyp;
	int smodel;
	double[] params;
	int rmodel;
	int alphaCats;
	double alpha;
	double fracInv;
	int numDataSets;
	int numSites;
	int treeTest;
	int numBootstraps;
	int format, distMethod, branchConstraint, dropSites, treeMethod;
	double[] freq;
	boolean mltree, mldist, evolve, disttree, rewrite;
	boolean userFreqs, optTree, optModel, useModel, siteLik;
	boolean weightedLS, dropGaps, jumble, bootstrap;
	int numTrees;  // specified when calling setMLTREE
	
	
	void setMLTREE(int numTrees)
	{
		mltree = true;
		mldist = false;
		evolve = false;
		disttree = false;
		rewrite = false;
		
		this.numTrees = numTrees;
	}

	void setMLDIST()
	{
		mltree = false;
		mldist = true;
		evolve = false;
		disttree = false;
		rewrite = false;
	}

	void setEVOLVE()
	{
		mltree = false;
		mldist = false;
		evolve = true;
		disttree = false;
		rewrite = false;
	}
	
	void setDISTTREE()
	{
		mltree = false;
		mldist = false;
		evolve = false;
		disttree = true;
		rewrite = false;
	}

	void setREWRITE()
	{
		mltree = false;
		mldist = false;
		evolve = false;
		disttree = false;
		rewrite = true;
	}

	void setOptions()
	{
		changeOptions();
	}
	
	
	//
	// Private stuff
	//
		
	private void resetParams()
	{
		for (int i = 0; i < params.length; i++)
		{
			params[i] = 1.0;
		}
	}	
	
	// Print "not possible" message
	private void notPossible()
	{
		System.out.println();
		System.out.println();
		System.out.println();
		System.out.println("This is not a possible option!");
	}

	// Print "unvalid input" message
	private void notValid()
	{
		System.out.println();
		System.out.println();
		System.out.println();
		System.out.println("This input is not valid!");
	}

	
	// Print current set of parameters
	private void printOptions()
	{
		System.out.println();
		System.out.println();
		System.out.println();
		System.out.println("GENERAL OPTIONS");
		if (rewrite)
		{
			System.out.print(" v               Sequence output format?  ");
			if (format == 0) System.out.println("PHYLIP Interleaved");
			if (format == 1) System.out.println("Plain (1 seq./line)");
			if (format == 2) System.out.println("CLUSTAL W");
			if (format == 3) System.out.println("PHYLIP sequential");
			System.out.print(" g                 Drop sites with gaps?  ");
			if (dropGaps) System.out.println("Yes");
			else System.out.println("No");

			System.out.print(" z       Constant/non-informative sites?  ");
			if (dropSites == 0) System.out.println("Keep sites");
			else if (dropSites == 1) System.out.println("Drop constant sites");
			else System.out.println("Drop non-informative sites");	
			
			System.out.print(" j                Jumble sequence order?  ");
			if (jumble) System.out.println("Yes");
			else System.out.println("No");
		
			
		}
		if (rewrite)
		{	
			System.out.print(" h                      Bootstrap sites?  ");
			if (bootstrap) System.out.println("Yes");
			else System.out.println("No");
		}
		
		
		if (mltree || evolve || mldist)
		{
			System.out.print(" d                   Sequence data type?  ");
			if (dtyp == DataType.NUCLEOTIDES) System.out.println("Nucleotides");
			if (dtyp == DataType.AMINOACIDS) System.out.println("Amino acids");
			if (dtyp == DataType.TWOSTATES) System.out.println("Binary states");
		}
		if ( (mltree || (mldist && useModel)) && ((rmodel == 1 || rmodel == 2) || (dtyp == DataType.NUCLEOTIDES && smodel < 4)))
		{
			System.out.print(" e            Estimate model parameters?  ");
			if (optModel) System.out.println("Yes");
			else System.out.println("No");
		}
		if (disttree)
		{
			System.out.print(" u                Reconstruction method?  ");
			if (distMethod == 0) System.out.println("Least-squares (on user trees)");
			if (distMethod == 1) System.out.println("Neighbor-Joining");
			if (distMethod == 2) System.out.println("UPGMA");
		}
		if (mltree || (disttree && distMethod == 0))
		{
			System.out.print(" o              Optimize branch lengths?  ");
			if (optTree) System.out.println("Yes");
		      	else System.out.println("No (use input branch lengths)");
			if (optTree)
			{
				System.out.print(" b           Branch lengths constraints?  ");
				if (branchConstraint == 0) System.out.println("None");
				if (branchConstraint == 1) System.out.println("Clock-like tree");
				if (branchConstraint == 2) System.out.println("Clock-like tree with dated tips");
			}
		}
		
		if (mltree)
		{	
			System.out.print(" i       Write out site log-likelihoods?  ");
			if (siteLik) System.out.println("Yes");
			else System.out.println("No");
		}
		
		if (mltree && numTrees > 1)
		{	
			System.out.print(" x                 Tree comparison test?  ");
			if (treeTest == 0) System.out.println("None");
			else if (treeTest == 1) System.out.println("Credible Sets");
			else if (treeTest == 2) System.out.println("Shimodaira-Hasegawa");
			else if (treeTest == 3) System.out.println("Kishino-Hasegawa");
			else if (treeTest == 4) System.out.println("All");
		}
		
		if (disttree && distMethod == 0)
		{
			System.out.print(" l                  Chi-square distance?  ");
			if (weightedLS) System.out.println("Weighted (Fitch-Margoliash 1967)");
			else System.out.println("Unweighted (Cavalli Sforza-Edward 1967)");
		}
		
		if (mldist)
		{
			System.out.print(" k     Compute observed or ML distances?  ");
			if (useModel) System.out.println("ML distances");
			if (!useModel) System.out.println("Observed distances");
		}
		

		if (evolve || rewrite)
		{
			System.out.println(" n           Number of output data sets?  " + numDataSets);
		}
		if (mldist || (mltree && optTree && numTrees == 1) )
		{
			System.out.println(" n                 Number of bootstraps?  " + numBootstraps);
		}
		
		if (mldist)
		{	
			System.out.print(" u           Tree reconstruction method?  ");
			if (treeMethod == 0) System.out.println("Neighbor-Joining");
			if (treeMethod == 1) System.out.println("UPGMA");
		}
		
		if (evolve)
		{
			System.out.println(" s                      Number of sites?  " + numSites);
		}
		if (mltree || evolve || (mldist && useModel))
		{
			System.out.println("SUBSTITUTION MODEL");
			System.out.print(" m                Model of substitution?  ");
			if (dtyp == DataType.NUCLEOTIDES)  // Nucleotides
			{
				if (smodel == NucleotideModelID.GTR) System.out.println("GTR (Lanave et al. 1984)");
				if (smodel == NucleotideModelID.TN) System.out.println("TN (Tamura-Nei 1993)");
				if (smodel == NucleotideModelID.HKY) System.out.println("HKY (Hasegawa et al. 1985)");
				if (smodel == NucleotideModelID.F84) System.out.println("F84 (Felsenstein 1984, PHYLIP)");
				if (smodel == NucleotideModelID.F81) System.out.println("F81 (Felsenstein 1981)");

				if (!optModel || evolve)
				{
					if (smodel == NucleotideModelID.GTR)
					{
						System.out.print(" r                       GTR parameters?  ");
						System.out.println(params[0] + ", " + params[1] + ", " +
							params[2] + ", " + params[3] + ", " + params[4]);			
					}							
					if (smodel == NucleotideModelID.TN || smodel == NucleotideModelID.HKY)
					{
						System.out.print(" t               Ts/Tv rate ratio kappa?  ");
						System.out.println(params[0]);
	
					}
					if (smodel == NucleotideModelID.F84)
					{
						System.out.print(" t               PHYLIP Ts/Tv parameter?  ");
						System.out.print(params[0]);
						System.out.println(" (expected Ts/Tv ratio)");
					}
					
					if (smodel == NucleotideModelID.TN)
					{
						System.out.print(" r            Y/R transition rate ratio?  ");
						System.out.println(params[1]);			
					}
				}			
			}
		
			if (dtyp == DataType.AMINOACIDS)  // Amino acids
			{
				if (smodel == AminoAcidModelID.DAYHOFF) System.out.println("Dayhoff (Dayhoff et al. 1978)");	
				if (smodel == AminoAcidModelID.JTT) System.out.println("JTT (Jones et al. 1992)");
				if (smodel == AminoAcidModelID.MTREV24) System.out.println("MTREV24 (Adachi-Hasegawa 1996)");
				if (smodel == AminoAcidModelID.BLOSUM62) System.out.println("BLOSUM62 (Henikoff-Henikoff 1992)");
				if (smodel == AminoAcidModelID.VT) System.out.println("VT (Mueller-Vingron 2000)");
				if (smodel == AminoAcidModelID.WAG) System.out.println("WAG (Whelan-Goldman 2000)");
				if (smodel == AminoAcidModelID.CPREV) System.out.println("CPREV (Adachi et al. 2000)");
			}
			if (dtyp == DataType.TWOSTATES) // Binary states
			{
				if (smodel == 0) System.out.println("F81 (Felsenstein 1981)");
			}
			
			if (dtyp == DataType.NUCLEOTIDES) System.out.print(" f               Nucleotide frequencies?  ");
			if (dtyp == DataType.AMINOACIDS) System.out.print(" f               Amino acid frequencies?  ");
			if (dtyp == DataType.TWOSTATES) System.out.print(" f                Two-state frequencies?  ");
			
			if (userFreqs)
			{
				System.out.println("Use specified values");
			}
			else if (evolve)
			{
				System.out.println("Not yet specified");
			}
			else
			{
				System.out.println("As observed in data set");
			}	
							
			System.out.print(" w          Model of rate heterogeneity?  ");
			if (rmodel == 0) System.out.println("Uniform rate");
			if (rmodel == 1) System.out.println("Gamma distributed rates");
			if (rmodel == 1)
			{
				if (!optModel)
				{
					System.out.print(" a   Gamma distribution parameter alpha?  ");
					if (alpha < 1.0)
						System.out.println(alpha + " (strong rate heterogeneity)");
					else
						System.out.println(alpha + " (weak rate heterogeneity)");
				}
				
				System.out.println(" c      Number of Gamma rate categories?  " + alphaCats);
			}
			if (rmodel == 2) System.out.println("Invariable sites model");
			if (rmodel == 2)
			{
				if (!optModel)
				{
					System.out.println(" a         Fraction of invariable sites?  " + fracInv);
				}
			}
		}
	
		System.out.println();
		System.out.print("Quit [q], confirm [y], or change [menu] settings: ");
	}

	// Change options interactively
	private void changeOptions()
	{	
		PushbackReader in = InputSource.openStdIn();
		FormattedInput fi = FormattedInput.getInstance();
		
		char c;
		do
		{
			printOptions();
			
			try
			{
				String str = fi.readLine(in, false);
				if (str.length() > 0)
				{
					c = str.charAt(0);
				}
				else
				{
					c = '\n';
				}
			}
			catch (Exception e)
			{
				System.out.println ("EXCEPTION");
				c = ' ';
			}
	
			switch (c)
			{
				case '\n':	break;
				
				case '\r':	break;
				
				case 'a':	if (useModel && !optModel && (rmodel == 1 || rmodel == 2))
						{
							if (rmodel == 1)
							{
								System.out.println();
								System.out.print("Gamma distribution parameter alpha:  ");
								try
								{
									alpha = fi.readDouble(in);
									fi.nextLine(in);
								}
								catch (Exception e)
								{
									alpha = 1.0;
								}
							}
							else
							{
								System.out.println();
								System.out.print("Fraction of invariable sites:  ");
								try
								{
									fracInv = fi.readDouble(in);
									fi.nextLine(in);
								}
								catch (Exception e)
								{
									fracInv = 0.0;
								}
								if (fracInv > 1.0) fracInv = 1.0;
								if (fracInv < 0.0) fracInv = 0.0;
							}
						}
						else
						{
							notPossible();
						}
						break;

				case 'b':	if ((mltree || (disttree && distMethod == 0)) && optTree)
						{
							branchConstraint++;
							if (branchConstraint > 2) branchConstraint = 0;
						}
						else
						{
							notPossible();
						}
						break;

				case 'c':	if (useModel && rmodel == 1)
						{
							System.out.println();
							System.out.print("Number of Gamma rate categories:  ");
							try
							{
								alphaCats = fi.readInt(in);
								fi.nextLine(in);
							}
							catch (Exception e)
							{
								alphaCats = 8;
							}
							
							if (alphaCats <= 1)
							{
								rmodel = 0;
							}
						}
						break;

				case 'd':	if (mltree || evolve || mldist)
						{
				
							dtyp = dtyp + 1;
							if (dtyp == 3) dtyp = 0;
						
							// Change model
							smodel = 0;
							freq = null;
							userFreqs = false;
							break;
						}
						else
						{
							notPossible();
						}
						break;

				case 'e':	if ((mltree || (mldist && useModel)) &&
							((rmodel == 1 || rmodel == 2) ||
							(dtyp == DataType.NUCLEOTIDES && smodel != NucleotideModelID.F81)))
						{
							if (optModel) optModel = false;
							else optModel = true;
						}
						else
						{
							notPossible();
						}
						break;

				
				case 'f':	if (useModel)
						{
							if (userFreqs)
							{
								userFreqs = false;
								freq = null;
							}
							else
							{
								//DataType dataType = DataTypeUtils.getInstance(dtyp);
								DataType dataType = DataType.Utils.getInstance(dtyp);
								
								System.out.println();
								System.out.print(dataType.getDescription());
								System.out.println(" frequencies (in %):");
								System.out.println();
							
								int numStates = dataType.getNumStates();
								freq = new double[numStates];
								userFreqs = true;
							
								try
								{
									double sum = 0.0;
									for (int i = 0; i < numStates-1; i++)
									{
										System.out.print("pi(" + 
											dataType.getChar(i) + ")=");
										freq[i] = fi.readDouble(in);
										fi.nextLine(in);
							
										freq[i] = freq[i]/100.0;
										sum = sum + freq[i];
										if (sum > 1.0 || freq[i] < 0.0)
											throw new Exception();
									}
									freq[numStates-1] = 1.0-sum;
								 
								}
								catch (Exception e)
								{
									userFreqs = false;
									freq = null;
									notValid();
								}
							}
									
						}
						else
						{
							notPossible();
						}
						break;

				case 'g':	if (rewrite)
						{
							if (dropGaps) dropGaps = false;
							else dropGaps = true;
						}
						else
						{
							notPossible();
						}
						break;

				case 'h':	if (rewrite)
						{
							if (bootstrap) bootstrap = false;
							else bootstrap = true;
						}
						else
						{
							notPossible();
						}
						break;

				case 'i':	if (mltree)
						{
							if (siteLik) siteLik = false;
							else siteLik = true;
						}
						else
						{
							notPossible();
						}
						break;

				case 'j':	if (rewrite)
						{
							if (jumble) jumble = false;
							else jumble = true;
						}
						else
						{
							notPossible();
						}
						break;

				
				case 'k':	if (mldist)
						{
							if (useModel) useModel = false;
							else useModel = true;
						}
						else
						{
							notPossible();
						}
						break;

				case 'l':	if (disttree && distMethod == 0)
						{
							if (weightedLS) weightedLS = false;
							else weightedLS = true;
						}
						else
						{
							notPossible();
						}
						break;


				case 'm':	if (mltree || evolve || (mldist && useModel))
						{
							smodel = smodel + 1;
							if (dtyp == DataType.NUCLEOTIDES)
							{
								if (smodel == NucleotideModelID.MODELCOUNT) smodel = 0; 
							}
							if (dtyp == DataType.AMINOACIDS)
							{
								if (smodel == AminoAcidModelID.MODELCOUNT) smodel = 0; 
							}
							if (dtyp == DataType.TWOSTATES)
							{
								if (smodel > 0) smodel = 0; 
							}
						}
						else
						{
							notPossible();
						}
						break;

				case 'n':	if (evolve || rewrite || mldist || (mltree && optTree && numTrees == 1))
						{
							System.out.println();
							
							if (mldist || (mltree && optTree && numTrees == 1))
							{
								System.out.print("Number of bootstraps:  ");
							
								try
								{
									numBootstraps = fi.readInt(in);
									fi.nextLine(in);
								}
								catch (Exception e)
								{
									numBootstraps = 1;
								}
							}
							else
							{
								System.out.print("Number of data sets:  ");
							
								try
								{
									numDataSets = fi.readInt(in);
									fi.nextLine(in);
								}
								catch (Exception e)
								{
									numDataSets = 1;
								}
		
							}
							
							
						}
						else
						{
							notPossible();
						}
						break;

						
				case 'o':	if (mltree || (disttree && distMethod == 0))
						{
							if (optTree) optTree = false;
							else optTree = true;
						}
						else
						{
							notPossible();
						}
						break;

				case 'q':	System.exit(1);
						break;

				case 'r':	if (useModel && dtyp == DataType.NUCLEOTIDES &&
							smodel == NucleotideModelID.TN && !optModel)
						{
							System.out.println();
							System.out.print("Y/R transition rate ratio:  ");
							try
							{
								params[1] = fi.readDouble(in);
								fi.nextLine(in);
							}
							catch (Exception e)
							{
								params[1] = 1.0;
							}
						}
						else if (useModel && dtyp == DataType.NUCLEOTIDES &&
							smodel == NucleotideModelID.GTR && !optModel)
						{
							try
							{
								System.out.println();
								System.out.print("REV parameter a:  ");
								params[0] = fi.readDouble(in);
								fi.nextLine(in);
								System.out.print("REV parameter b:  ");
								params[1] = fi.readDouble(in);
								fi.nextLine(in);
								System.out.print("REV parameter c:  ");
								params[2] = fi.readDouble(in);
								fi.nextLine(in);
								System.out.print("REV parameter d:  ");
								params[3] = fi.readDouble(in);
								fi.nextLine(in);
								System.out.print("REV parameter e:  ");
								params[4] = fi.readDouble(in);
								fi.nextLine(in);
							}
							catch (Exception e)
							{
								resetParams();
							}
						}
						else
						{
							notPossible();
						}
						break;

				case 's':	if (evolve)
						{
							System.out.println();
							System.out.print("Number of sites:  ");
							try
							{
								numSites = fi.readInt(in);
								fi.nextLine(in);
							}
							catch (Exception e)
							{
								numSites = 100;
							}
						}
						else
						{
							notPossible();
						}
						break;

				case 't':	if (useModel && dtyp == DataType.NUCLEOTIDES &&
						(smodel == 1 || smodel == 2 || smodel == 3) && !optModel)
						{
							System.out.println();
							
							if (smodel != 3)
								System.out.print("Transition/transversion rate ratio:  ");
							else
								System.out.print("PHYLIP Ts/Tv parameter:  ");
							try
							{
								params[0] = fi.readDouble(in);
								fi.nextLine(in);
							}
							catch (Exception e)
							{
								params[0] = 1.0;
							}
						}
						else
						{
							notPossible();
						}
						break;

				case 'u':	if (disttree)
						{
							distMethod++;
							if (distMethod > 2) distMethod = 0;
						}
						else if (mldist)
						{
							treeMethod++;
							if (treeMethod > 1) treeMethod = 0;
						}
						else
						{
							notPossible();
						}
						break;


				case 'v':	if (rewrite)
						{
							format++;
							if (format > 3) format = 0;
						}
						else
						{
							notPossible();
						}
						break;

				case 'w':	if (useModel)
						{
							rmodel = rmodel + 1;
							if (rmodel > 2) rmodel = 0;
						}
						else
						{
							notPossible();
						}
						break;

				case 'x':	if (mltree && numTrees > 1)
						{
							treeTest++;
							if (treeTest > 4) treeTest = 0;
						}
						else
						{
							notPossible();
						}
						break;

				case 'y':	break;
				
				case 'z':	if (rewrite)
						{
							dropSites++;
							if (dropSites > 2) dropSites = 0;
						}
						else
						{
							notPossible();
						}
						break;
			
						
				default:	notPossible();
						break;
			}
			
				
		} while (c != 'y');
	}
 }

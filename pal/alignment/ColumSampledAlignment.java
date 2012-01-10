package pal.alignment;

import pal.math.UrnModel;

public class ColumSampledAlignment  extends AbstractAlignment{
	
	//
	// Public stuff
	//

	/**
	 * Constructor
	 *
	 * @param raw original alignment
	 */
	public ColumSampledAlignment(Alignment raw, int[] idx)
	{
		rawAlignment = raw;

		numSeqs = raw.getSequenceCount();
		idGroup = raw;
		numSites = idx.length;
		setDataType(raw.getDataType());
		this.alias=idx;
	}

	// Implementation of abstract Alignment method

	/** sequence alignment at (sequence, site) */
	public char getData(int seq, int site)
	{
		return rawAlignment.getData(seq, alias[site]);
	}



	//
	// Private stuff
	//


	private Alignment rawAlignment;
	private int[] alias;

}

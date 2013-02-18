package app;


import java.util.Map;
import java.util.HashMap;
import java.util.TreeMap;
import java.util.Map.Entry;
import java.util.Iterator;
import java.util.List;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.Collections;
import java.util.Random;
import java.util.Set;
import java.util.HashSet;
import java.util.TreeSet;
import java.util.Arrays;
import java.text.DecimalFormat;


import segment.SegmentList;
import segment.Segment;


public class AnnotationData {

	private double fdrmaxp = 0;                   // largest p value considered in the FDR summary
	private int verbose = 0;
	private int bucketsize = 1;

	static final String sampleKeyString = "##Samples";
	static final String observationKeyString = "##Observed";
	static final String annotationKeyString = "##Annotation";
	static final String segmentKeyString = "##Segments";

	public AnnotationData() {

		workspace = new HashMap();
		segments = new HashMap();
		ids = new HashMap();
		annotation = new HashMap();
		synonyms = new HashMap();
		original = new HashMap();
		subsets = new HashMap();

		// Load class file, to make us independent of jar file changes
		FalseDiscoveryStats fds = new FalseDiscoveryStats();
		fds.numOverrepresented = 0;  // to get rid of 'unused variable' warning

	}


	void intersectWorkspaces( int verbose ) {
		// This intersects workspaces when one or more are defined on the 'original', and one or more is
		// defined on the isochores

		Set originals = new HashSet(0);

		Iterator iter = workspace.keySet().iterator();
		while (iter.hasNext()) {
			String chrom = (String)iter.next();
			// is it an isochore identifier?
			if (original.containsKey( chrom )) {
				String orig = (String)original.get( chrom );
				// double check that original is different from current
				if (!orig.equals(chrom)) {
					// does original exist as workspace?
					SegmentList origworkspace;
					if (!workspace.containsKey( orig )) {
						System.out.println("# Warning: isochore "+chrom+" lives on chromosome "+orig+" for which no workspace is defined; assuming empty.");
						origworkspace = new SegmentList();
					} else {
						origworkspace = (SegmentList)workspace.get( orig );
						// Store original, to remove it later
						originals.add( orig );
					}
					// intersect this workspace with the original
					((SegmentList)workspace.get(chrom)).intersect( origworkspace );
					if (verbose > 1) {
						System.out.println("# Intersecting workspace "+chrom+" with "+orig);
					}
				}
			}
		}

		// Remove workspaces
		iter = originals.iterator();
		while (iter.hasNext()) {
			String chrom = (String)iter.next();
			workspace.remove( chrom );
			if (verbose > 1) {
				System.out.println("# Removing now-redundant workspace "+chrom);
			}
		}
	}


	void checkMutualWorkspaceIntersection() {

		// Checks that isochores are mutually non-intersecting

		Iterator iter = workspace.keySet().iterator();
		while (iter.hasNext()) {
			String iso = (String)iter.next();
			String orig = (String)original.get( iso );
			Iterator iter2 = workspace.keySet().iterator();
			while (iter2.hasNext()) {
				String iso2 = (String)iter2.next();
				if (!iso2.equals(iso)) {
					String orig2 = (String)original.get( iso2 );
					if (orig2.equals(orig)) {
						long intersection = ((SegmentList)workspace.get(iso)).intersection( (SegmentList)workspace.get(iso2) );
						if (intersection > 0) {
							System.out.println("# Warning -- isochores "+iso+" and "+iso2+" intersect by "+intersection+" nucleotides.");
						}
					}
				}
			}
		}
	}



	void addDefaultSynonyms( int verbose ) {

		Iterator iter = workspace.keySet().iterator();
		boolean first = true;
		while (iter.hasNext()) {
			String chrom = (String)iter.next();
			if (!original.containsKey( chrom ) ) {
				if (first) {
					if (verbose > 1) {
						System.out.print("# Adding default synonyms for: ");
						first = false;
					}
				}
				List list = new ArrayList(1);
				list.add( chrom.intern() );
				synonyms.put( chrom, list );
				original.put( chrom, chrom );
				if (verbose > 1) {
					System.out.print(" "+chrom);
				}
			}
		}
		if (!first && (verbose>1)) {
			System.out.println("");
		}
	}


	void checkSubsets( int verbose ) {
		// Adds default subset if none were defined, and checks that all subset ids map to existing synonyms

		if (subsets.isEmpty()) {
			List ss = new FastArrayList(0);
			ss.addAll( original.keySet() );
			subsets.put( new String(""), ss );
			if (verbose > 0) {
				System.out.println("# Adding all synonyms as default subset");
			}
		}

		Iterator iter = subsets.keySet().iterator();
		while (iter.hasNext()) {
			String subsetid = (String)iter.next();
			Iterator syniter = ((List)subsets.get(subsetid)).iterator();
			while (syniter.hasNext()) {
				String syn = (String)syniter.next();
				if (!original.containsKey( syn )) {
					System.out.println("Error: Subset '"+subsetid+"' refers to undefined synonym '"+syn+"'");
				}
			}
		}    	
	}



	void checkChroms() {

		Iterator iter = ids.entrySet().iterator();
		while (iter.hasNext()) {
			Entry e = (Entry)iter.next();
			Location loc = (Location)e.getValue();
			List syns = (List)synonyms.get( loc.chrom );
			if (syns == null) {
				System.out.println("# Warning: ID "+e.getKey()+" is located on chromosome "+loc.chrom+" for which no workspace is defined.\n#  Adding empty workspace.");
				workspace.put( loc.chrom, new SegmentList() );
				addDefaultSynonyms( 2 );
			}
		}
		iter = segments.keySet().iterator();
		while (iter.hasNext()) {
			String chrom = (String)iter.next();
			List syns = (List)synonyms.get( chrom );
			if (syns == null) {
				System.out.println("# Warning: Segments found for chromosome "+chrom+" for which no workspace is defined.\n#  Adding empty workspace.");
				workspace.put( chrom, new SegmentList() );
				addDefaultSynonyms( 2 );
			}
		}
	}



	int undefinedIDs( int verbose ) {
		// Build a set of undefined IDs
		Set undefIDs = new HashSet();
		Iterator annIter = annotation.keySet().iterator();
		while (annIter.hasNext()) {
			long[] idlist = (long[])annotation.get( (String)annIter.next() );
			for (int i=0; i<idlist.length; i++) {
				Long id = new Long( idlist[i] );
				if (!ids.keySet().contains( id )) {
					undefIDs.add( id );
				}
			}
		}
		if (verbose > 2) {
			if (undefIDs.size() > 0) {
				System.out.println("# Found undefined IDs (ignoring):\n#");
				Iterator iter = undefIDs.iterator();
				while (iter.hasNext()) {
					Long id = (Long)iter.next();
					System.out.print(id + "\t");
				}
				System.out.println("");
			}
		}
		return undefIDs.size();
	}



	int purgeUnusedIDs() {
		// Build a set of unused IDs
		Set unusedIDs = new HashSet( ids.keySet() );
		Iterator annIter = annotation.keySet().iterator();
		while (annIter.hasNext()) {
			long[] ids = (long[])annotation.get( (String)annIter.next() );
			for (int i=0; i<ids.length; i++) {
				Long id = new Long( ids[i] );
				// Remove the used ID - if it exists
				unusedIDs.remove( id );
			}
		}

		// Loop over IDs and remove unused ones
		Iterator idIter = unusedIDs.iterator();
		while (idIter.hasNext()) {
			Long id = (Long)idIter.next();
			ids.remove( id );
		}

		return unusedIDs.size();

	}



	int removeIdsOutsideWorkspace() {

		Iterator iter = ids.entrySet().iterator();
		Set outside = new HashSet();
		boolean firstIteration = true;
		while (iter.hasNext()) {
			Entry e = (Entry)iter.next();
			Long id = (Long)e.getKey();
			Location loc = (Location)e.getValue();
			long intersection = 0;
			List syns = (List)synonyms.get( loc.chrom );
			for (int idx=0; idx<syns.size(); idx++) {
				String isoch = (String)syns.get(idx);
				SegmentList segList = (SegmentList)workspace.get( isoch );
				if (segList == null) {
					if (firstIteration) {
						System.out.println("Warning - no segments found for workspace '"+isoch+"'");
					}
				} else {
					intersection += segList.intersection( loc.getSegments() );
				}
			}
			if (intersection == 0) {
				outside.add( id );
			}
			firstIteration = false;
		}
		// Loop over IDs and remove unused ones from id list
		Iterator idIter = outside.iterator();
		while (idIter.hasNext()) {
			Long id = (Long)idIter.next();
			ids.remove( id );
		}
		// Return the number of ids removed
		return outside.size();
	}



	void putIds( Long newId, Location newLoc ) {

		ids.put( newId, newLoc );

	}



	class IdArrayComparator implements Comparator {

		public int compare( Object a, Object b) {

			long[] p = (long[])a;
			long[] q = (long[])b;
			if (p.length != q.length) {
				return p.length - q.length;
			}
			for (int i=0; i<p.length; i++) {
				if (p[i] != q[i]) {
					if (p[i]>q[i]) {
						return 1;
					} else {
						return -1;
					}
				}
			}
			return 0;    		
		}

	}


	int removeEmptyAnnotations( int verbose ) {

		int removed = 0;
		Iterator annIter = new ArrayList( annotation.entrySet() ).iterator();
		while (annIter.hasNext()) {
			Entry e = (Entry)annIter.next();
			long[] idarray = (long[])e.getValue();
			// remove id values that are no (longer) existent
			int ptr = 0;
			for (int i=0; i<idarray.length; i++) {
				Long id = new Long( idarray[i] );
				if (ids.containsKey( id )) {
					// Keep this one
					idarray[ptr] = idarray[i];
					ptr += 1;
				}
			}
			if (ptr == 0) {
				if (verbose > 1) {
					if ((removed < 100) || (verbose > 2)) {
						System.out.println("# Removing empty annotation: "+(String)e.getKey());
					}
					if ((removed == 100) && (verbose == 2)) {
						System.out.println("#  (More annotations removed; further output suppressed.)");
					}
				}
				annotation.remove( e.getKey() );
				removed += 1;
			} else {
				if (ptr != idarray.length) {
					long[] newids = new long[ptr];
					System.arraycopy( idarray, 0, newids, 0, ptr );
					annotation.put( e.getKey(), newids );
				}
			}
		}

		return removed;
	}


	int removeDuplicatedAnnotations( int verbose ) {

		Iterator annIter = annotation.values().iterator();
		while (annIter.hasNext()) {
			Arrays.sort( (long[])annIter.next() );
		}
		List idArrays = new ArrayList( annotation.values() );
		IdArrayComparator idArrayComparator = new IdArrayComparator();
		Collections.sort( idArrays, idArrayComparator );
		boolean removing = false;
		int removed = 0;
		for (int i=1; i<idArrays.size(); i++) {
			long[] previous = (long[])idArrays.get(i-1);
			long[] current = (long[])idArrays.get(i);
			if (Arrays.equals( previous, current )) {
				// Now find 'm again...
				String keyPrev = null, keyCurr = null;
				annIter = annotation.entrySet().iterator();
				while (annIter.hasNext()) {
					Entry e = (Entry)annIter.next();
					if (current == (long[])e.getValue()) {
						keyCurr = (String)e.getKey();
					}
					if (previous == (long[])e.getValue()) {
						keyPrev = (String)e.getKey();
					}
				}
				if (!removing) {
					removing = true;
					if (verbose > 1) {
						System.out.println("# Found duplicated annotation: "+keyPrev);
					}
				}
				if (verbose > 1) {
					System.out.println("#   Removing:                  "+keyCurr);
				}
				annotation.remove( keyCurr );
				removed += 1;
			} else {
				removing = false;
			}
		}
		return removed;
	}


	public void reportSegmentIntersection( int verbose ) {

		if (verbose == 0) {
			return;
		}

		// total intersection of segments with base chromosome
		HashMap chromIntersection = new HashMap(0);
		int totalSegs = 0;     // total number of segments
		int numSegs = 0;       // total number of segments in workspace
		int numFulls = 0;      // number of segments fully within workspace
		boolean warnOverlap;

		// loop over all workspaces
		Iterator iter = workspace.entrySet().iterator();
		while (iter.hasNext()) {
			// get workspace segments
			Entry e = (Entry)iter.next();
			SegmentList workList = (SegmentList)e.getValue();
			// get chromosome (=isochore) and base chromosome
			String chrom = (String)e.getKey();
			String basechrom = (String)original.get( chrom );
			if (basechrom == null) {
				basechrom = chrom;
			}
			// enter map for intersection data
			HashMap segdata;
			if (!chromIntersection.containsKey( basechrom) ) {
				segdata = new HashMap(0);
				chromIntersection.put( basechrom, segdata );
				// first time we see this chromosome - accumulate segment total
				totalSegs += ((SegmentList)segments.get(basechrom)).numSegments();
			} else {
				segdata = (HashMap)chromIntersection.get( basechrom );
			}
			// get input segment list
			SegmentList segList = (SegmentList)segments.get(basechrom);   // segment lists are defined on base chrom
			if (segList == null) {
				segList = new SegmentList();
			}
			// intersection with workspace
			List segs = segList.intersectingSegments( workList );
			Iterator segIter = segs.iterator();
			warnOverlap = true;
			long isochoreOverlap = 0;
			while (segIter.hasNext()) {                     
				Segment seg = (Segment)segIter.next();
				long segLength = seg.length();
				long intersectLength = workList.intersection( seg );
				long totLength = 0;
				if (segdata.containsKey( seg )) {
					totLength = ((Integer)segdata.get( seg )).intValue();
				}
				if ((totLength == 0) && (intersectLength > 0)) {
					// first-time hit - add to total count
					numSegs += 1;
				}
				// we assume that the isochores are disjunct; the program gives a warning if this is not so
				totLength += intersectLength;
				segdata.put( seg, new Integer( (int)totLength ) );
				if ((totLength >= segLength)&&(intersectLength > 0)) {
					numFulls += 1;
				};
				if (totLength > segLength && warnOverlap) {
					System.out.println("# Warning - total overlap of isochores with segment is longer than segment in chrom "+chrom+" (isochore problem?)  Segment counts unreliable.");
					System.out.println("#  (Offending segment: "+seg+" )");
					// do not generate any further warnings for this isochore
					warnOverlap = false;
				};
				isochoreOverlap += intersectLength;
			}
			if (verbose > 3) {
				System.out.println("# ReportIntersection: Segment overlap with workspace in isochore "+chrom+": "+isochoreOverlap);
			}
		}

		System.out.println("# Number of segments: "+totalSegs);
		System.out.println("# Number of segments within workspace: "+numSegs);
		System.out.println("# Number of segments fully within workspace: "+numFulls);
		System.out.println("# (Note: Segments that touch or overlap are merged before counting),");
	}


	/**
	 * Assign input segments to current seg;
	 * calculate intersection with workspace;
	 * calculate length distribution
	 */
	public void initialize(Random random, double baseline, int numIterations, double fdrmaxp, int verbose, 
			int histsize, int bucketsize, boolean keepOverhang, boolean dumpSamples, boolean intersected) {

		this.verbose = verbose;
		this.fdrmaxp = fdrmaxp;
		this.bucketsize = bucketsize;

		// Assign random number generator to private field
		this.random = random;

		// Assign observed set IGS to current set
		currentSegments = segments;

		workIgsLength = new HashMap(0);		
		largeWorkspace = new HashMap(0);
		lengthDistr = new HashMap(0);
		currentWorkSegs = new HashMap(0);
		totalLength = new HashMap(0);
		annotationStats = new TreeMap();	// better iteration and memory behaviour
		idIntersection = new HashMap(0);
		segmentListSamplers = new HashMap(0);
		fullWorkspaces = new HashMap(0);

		// Reset total IGS intersection for all subsets
		Iterator subsetiter = subsets.keySet().iterator();
		while (subsetiter.hasNext()) {
			String subset = (String)subsetiter.next();
			workIgsLength.put(subset, new Long(0));
		}

		// Collect statistics of over-sized segments
		int numSegs = 0;
		int numOversize = 0;
		long totOversize = 0;
		long maxOversize = 0;

		// Loop over all workspace isochores
		Iterator iter = workspace.entrySet().iterator();
		while (iter.hasNext()) {
			Entry e = (Entry)iter.next();
			SegmentList workList = (SegmentList)e.getValue();
			String chrom = (String)e.getKey();
			String basechrom = (String)original.get( chrom );
			if (basechrom == null) {
				System.out.println("# Warning: no chromosome found for isochore "+chrom+" referred to in workspace.  Syn file problem?  Proceeding...");
				basechrom = chrom;
			}
			if (!intersected) {
			    // build up full work list
			    SegmentList fullWorkspace = (SegmentList)fullWorkspaces.get( basechrom );
			    if (fullWorkspace == null) fullWorkspace = new SegmentList();
			    fullWorkspace.union( workList );
			    fullWorkspaces.put( basechrom, fullWorkspace );
			}
			SegmentList segList = (SegmentList)currentSegments.get(basechrom);   //!
			if (segList == null) {
				System.out.println("# Warning: No segments for chromosome "+basechrom+" found.  Assuming empty.");
				segList = new SegmentList();
				currentSegments.put( basechrom, segList );
			}
			currentSegments.put( chrom, segList );  //!

			// Compute length distribution of IGS intersecting workspace isochore
			int[] histogram = new int[histsize];
			lengthDistr.put( chrom, histogram );
			List segs = segList.intersectingSegments( workList );
			Iterator segIter = segs.iterator();
			while (segIter.hasNext()) {                     
				Segment seg = (Segment)segIter.next();
				long length = seg.length();
				if (!keepOverhang) {
					length = segList.intersection( seg );
				}
				numSegs += 1;
				if (length+bucketsize-1 < histsize * bucketsize) {
					histogram[(int)((length+bucketsize-1)/bucketsize)] += 1;
				} else {
				    	histogram[histsize-1] += 1;
					numOversize += 1;
					totOversize += length;
					if (length > maxOversize) {
						maxOversize = length;
					}
					// Too important to ignore, sometimes
					throw new Error("Error: Large segment found, doesn't fit in histogram.  Increase histogram or bucket size.");
				}
			}

			// Get mean and variance of segment length distribution
			double mean = histogramMean( histogram ) * bucketsize;
			double var = histogramVar( histogram ) * bucketsize * bucketsize;
			int extension = (int)( mean + 4.0*Math.sqrt(var) - 0.5 );

			// Get enlarged workspace
			SegmentList largews = new SegmentList( workList );
			largews.extend( extension, 0 );
			largeWorkspace.put( chrom, largews );
			SegmentListSampler sls = new SegmentListSampler( largews, random );
			segmentListSamplers.put( chrom, sls );

			// Initialize currentWorkSeg (intersection with workspace)
			SegmentList currentWorkSeg = new SegmentList( workList );
			currentWorkSeg.intersect( segList );
			currentWorkSegs.put( chrom, currentWorkSeg );

			// Update workIgsLength
			subsetiter = subsets.keySet().iterator();
			while (subsetiter.hasNext()) {
				String subset = (String)subsetiter.next();
				List syns = (List)subsets.get( subset );
				if (syns.contains( chrom )) {
					// Current isochore 'chrom' is contained within 'subset'
					long intersection = ((Long)workIgsLength.get( subset )).longValue() + currentWorkSeg.length();
					workIgsLength.put( subset, new Long(intersection) );
				}
			}
			
			// Calculate total length
			long totlen = currentWorkSeg.length();
			totalLength.put( chrom, new Long(totlen) );

			if (verbose>3) {
				System.out.println("# Initialize: Segment overlap with workspace in isochore "+chrom+": "+totlen);
			}
			if (verbose > 4) {
				currentWorkSeg.dump("##SegWork\t"+chrom);
			}

		}

		if (numOversize > 0) {
			System.out.println("# Warning: Found "+numOversize+" segments over "+histsize*bucketsize+" nt (of "+numSegs+"; max "+maxOversize+", mean "+totOversize/numOversize+")");
		}

		// Loop over all annotations; initialize annotationStats
		int keep = (int)(numIterations * fdrmaxp);
		Statistics.setKeep( keep );
		if (verbose>0) {
			System.out.println("# Keeping "+keep+" samples (for "+annotation.entrySet().size()+" annotations) for FDR statistics");
		}
		/*
		if (dumpSamples) {
		    Statistics.keepAll( numIterations );
		    if (verbose>0) {
			System.out.println("# Keeping a further "+numIterations+" samples for dumpSamples");
		    }
		}
		 */
		if (dumpSamples) {
			System.out.print(annotationKeyString);
		}
		// combine subset and annotations; also print header line for the sample dumps
		Iterator ssiter = subsets.keySet().iterator();
		while (ssiter.hasNext()) {
			String subset = (String)ssiter.next();
			iter = annotation.entrySet().iterator();
			while (iter.hasNext()) {
				String ann = (String)(((Entry)iter.next()).getKey());
				String annid = makeAnnId( subset, ann );
				if (dumpSamples) {
					System.out.print("\t" + annid);
				}
				annotationStats.put( annid, new Statistics( annid ) );
			}
		}
		if (dumpSamples) {
			System.out.println();
		}
		Statistics.setBaseline( baseline );

	}


	String makeAnnId( String subset, String ann ) {

		if (subset == "") {
			return ann;
		} else {
			return subset + "." + ann;
		}

	}


	void calculateIdintersection(boolean printWarning, String subset) {

		// Calculate idIntersection
		List subsetsyns = (List)subsets.get( subset );
		Iterator iter = ids.entrySet().iterator();
		Set notFound = new HashSet();
		while (iter.hasNext()) {
			Entry e = (Entry)iter.next();
			Long id = (Long)e.getKey();
			Location loc = (Location)e.getValue();
			long intersection = 0;
			List syns = (List)synonyms.get( loc.chrom );   //!
			for (int idx=0; idx<syns.size(); idx++) {
				String chrom = (String)syns.get(idx);
				if (subsetsyns.contains( chrom )) {
					// This synonym is part of the subset, so include it
					SegmentList segList = (SegmentList)currentWorkSegs.get( chrom );
					if (segList == null) {
						if (!notFound.contains(chrom)) {
							if (printWarning) {
								System.out.println("Warning - no segments found for '"+chrom+"'");
							}
							notFound.add(chrom);
						}
					} else {
						intersection += segList.intersection( loc.getSegments() );
					}
				}
			}
			idIntersection.put( id, new Long(intersection) );
		}
	}


	void calculate(int iteration, boolean dumpSamples, String subset, boolean externallysampled) {

		DecimalFormat df = new DecimalFormat();
		if (dumpSamples) {
			df.applyPattern("#.0000000");
			if (iteration < 0) {
				System.out.print(observationKeyString);
			} else {
				System.out.print(sampleKeyString);
			}
		}
		// Loop over all annotations, and update statistics
		Iterator iter = annotation.entrySet().iterator();
		while (iter.hasNext()) {
			Entry e = (Entry)iter.next();
			String ann = (String)e.getKey();
			String annid = makeAnnId( subset, ann );
			long[] ids = (long[])e.getValue();
			// Calculate number of nucs within segments intersecting with workspace, annotated as ann and within subset.
			long total = 0;
			for (int i=0; i<ids.length; i++) {
				Long intersect = (Long)idIntersection.get( new Long(ids[i]) );
				if (intersect != null) {
					total += intersect.longValue();
				}
			}
			// Update statistics
			double sample;
			if (externallysampled) {
			  //System.out.printf("externally:\t%d\t%d\t%d\t\n", total, ids.length, (Long)sampleIgsLength.get(subset));
			    sample = (0.0+total)/((Long)sampleIgsLength.get(subset)).longValue();
			} else {
			    sample = (0.0+total)/((Long)workIgsLength.get(subset)).longValue();
			    System.out.printf("## %s\t%d\t%d\t%d\t%f\n", annid, total, ids.length,
					      (Long)workIgsLength.get(subset), sample );
			}
			((Statistics)annotationStats.get( annid )).addSample( sample, iteration );
			if (dumpSamples) {
				System.out.print( "\t" + df.format(sample) );
			}
		}
	}



	void randomSample(boolean dumpSegments) {

		// Remove previous sample
		currentSegments = new HashMap(0);

		// Loop over all chromosomes (i.e. isochores)
		// Note: this assigns segments to isochores, not chromosomes as in input file!
		Iterator iter = workspace.entrySet().iterator();
		while (iter.hasNext()) {
			Entry e = (Entry)iter.next();
			String chrom = (String)e.getKey();
			SegmentList ws = (SegmentList)e.getValue();
			SegmentList segList = new SegmentList();
			HistogramSampler hs = new HistogramSampler( (int[])lengthDistr.get(chrom), bucketsize );
			SegmentListSampler sls = (SegmentListSampler)segmentListSamplers.get( chrom );

			long remaining = ((Long)totalLength.get(chrom)).longValue();
			long trueRemaining = remaining;	
			List newsegs = new ArrayList(1000);
			SegmentList intersectedSegments = new SegmentList();       // contains intersection of sampled segments with workspace

			int unsuccessfulRounds = 0;
			int maxUnsuccessfulRounds = 20;

			while (trueRemaining > 0) {
			    
			        // Sample a segment length from the histogram
			        int len = hs.sample(random);
			        if (len == 0) {
				    throw new Error("Sampled 0 from length histogram!");
				}

				// If we can now potentially get required number of nucleotides, recompute
				// the required amount since we might need more due to overlap
				if (remaining <= len) {
				    segList.addList( newsegs );
				    newsegs.clear();
				    intersectedSegments = new SegmentList( ws );
				    intersectedSegments.intersect( segList );
				    remaining = ((Long)totalLength.get(chrom)).longValue() - intersectedSegments.length();
				    if (trueRemaining == remaining) {
					unsuccessfulRounds += 1;
				    } else {
					trueRemaining = remaining;
				    }
				}

				// deal with overshoot
				if (trueRemaining < 0) {
				    long p = new SegmentListSampler( segList, sls.getRandom() ).sample();
				    segList.trim( p, -trueRemaining );
				}

				// sample a position until we get a nonzero overlap (but don't bother if len==0)
				long p = sls.sample();
				Segment seg = new Segment(p, p+len);
				long intersect = ws.intersection( seg );
				while (intersect == 0 && len > 0) {
				    p = sls.sample();
				    seg.left = p;
				    seg.right = p+len;
				    intersect = ws.intersection( seg );
				};
				
				// If intersect > remaining we could well add too much (but we're not sure, since
				// there may be overlap with previously sampled segments).  Make sure we don't
				// overshoot the mark.  However, if we can't increase the amount of overlap in a
				// reasonable amount of time, don't bother
				if (remaining < intersect && unsuccessfulRounds < maxUnsuccessfulRounds) {
				    if (ws.intersection( new Segment(p,p+1) ) == 1) {
					// Left end overlaps with workspace
					seg.right = seg.left + remaining;
				    } else {
					// Right end does?
					seg.left = seg.right - remaining;
				    }
				    intersect = ws.intersection( seg );
				}

				if (trueRemaining > 0) {
				    newsegs.add( seg );
				    remaining -= intersect;
				}

			}

			// Store this sample
			// currentSegments.put(chrom, segList);                 // no need to keep these
			currentSegments.put(chrom, null);                       // reduce memory usage
			currentWorkSegs.put(chrom, intersectedSegments);

			if (dumpSegments) {
			    intersectedSegments.dump(segmentKeyString + "\t" + chrom);
			}

		}
	}


    // Alternative sampler.  This one does not break up sampled segments, but may not sample the exact number
    // of bases in each GC category.  Therefore, "subset" results may be less reliable.  Whole-chromosome results
    // are more reliable than those using the other sampler, particularly when the isochores are defined on scales
    // that are similar or smaller than the segments themselves, since the previous sampler would show reduced
    // variance in this case.

    int isosampler( List lengths ) {

	long total = 0;
	int idx = 0;
	for (; idx<lengths.size(); idx++) {
	    total += (Long)lengths.get(idx);
	}
	if (total==0) return -1;
	long sample = Math.abs(random.nextLong()) % total;
	for (idx=0;; idx++) {
	    sample -= (Long)lengths.get(idx);
	    if (sample<0) return idx;
	}
    }


    void randomSample2(boolean dumpSegments) {
    
		// Remove previous sample
		currentSegments = new HashMap(0);

		// Loop over all chromosomes (NOT isochores)
		Iterator iter = synonyms.entrySet().iterator();
		while (iter.hasNext()) {
		    Entry e = (Entry)iter.next();
		    String chromosome = (String)e.getKey();
		    SegmentList fullworklist = (SegmentList)fullWorkspaces.get(chromosome);
		    List isochores = (List)e.getValue();
		    // calculate total required number of nucleotides, and required per isochore
		    List isoremaining = new ArrayList(isochores);
		    List isohistsamplers = new ArrayList(isochores);
		    long remaining = 0;
		    int idx = 0;
		    for (; idx < isochores.size(); idx++) {
			String isochore = (String)isochores.get(idx);
			Long isolen = (Long)totalLength.get( isochore );
			isoremaining.set( idx, isolen );
			remaining += isolen;
			HistogramSampler hs = new HistogramSampler( (int[])lengthDistr.get(isochore), bucketsize );
			isohistsamplers.set( idx, hs );
		    }
		    // initialization
		    List isototal = new ArrayList( isoremaining );          // to reinitialize the isoremaining array
		    SegmentList segList = new SegmentList();		    // sampled segments, not intersected
		    SegmentList intersectedSegments = new SegmentList();    // sampled segments intersected with workspace
		    long trueremaining = remaining;
		    long totalrequired = remaining;
		    int unsuccessfulRounds = 0;                             // to bail out, in case of densly packed chromosomes
		    int maxUnsuccessfulRounds = 20;
		    List newsegs = new ArrayList(1000);
		    // start loop proper
		    while (trueremaining > 0) {
		   
			// sample an isochore
			int iso = isosampler( isoremaining );
			if (iso == -1) {
			    // all counts went down to 0 -- this happens when there is overlap among
			    // the sampled segments.  Simply reset the counters to the original values
			    // to keep sampling at the right proportions
			    isoremaining.clear();
			    isoremaining.addAll( isototal );
			    iso = isosampler( isoremaining );
			    if (iso == -1) iso=0;  // failsafe
			}

			int len = 0;
			
			// sample a segment length
			len = ((HistogramSampler)isohistsamplers.get( iso )).sample( random );
			if (len == 0) {
			    throw new Error("Sampled 0 from length histogram!");
			}

			// If we can now potentially get required number of nucleotides, recompute
			// the required amount since we might need more due to overlap of sampled segments
			if (remaining <= len) {
			    segList.addList( newsegs );
			    newsegs.clear();
			    intersectedSegments = new SegmentList( fullworklist );
			    intersectedSegments.intersect( segList );
			    remaining = totalrequired - intersectedSegments.length();
			    if (trueremaining == remaining) {
				unsuccessfulRounds += 1;
			    } else {
				trueremaining = remaining;
			    }
			}

			// deal with overshoot
			if (trueremaining < 0) {
			    long position = new SegmentListSampler( segList, random ).sample();
			    segList.trim( position, -trueremaining );
			}

			// sample a position until we get a nonzero overlap (but don't bother if len==0)
			SegmentListSampler sls = (SegmentListSampler)segmentListSamplers.get( isochores.get(iso) );
			long p = sls.sample() - (len/2);
			Segment seg = new Segment(p, p+len);
			long intersect = fullworklist.intersection( seg );
			while (intersect == 0 && len>0) {
			    p = sls.sample();
			    seg.left = p;
			    seg.right = p+len;
			    intersect = fullworklist.intersection( seg );
			}
			// If intersect > remaining we could well add too much (but we're not sure, since
			// there may be overlap with previously sampled segments).  Make sure we don't
			// overshoot the mark.  However, if we can't increase the amount of overlap in a
			// reasonable amount of time, don't bother
			if (remaining < intersect && unsuccessfulRounds < maxUnsuccessfulRounds) {
			    if (fullworklist.intersection( new Segment(p,p+1) ) == 1) {
				// Left end overlaps with workspace
				seg.right = seg.left + remaining;
			    } else {
				// Right end does?
				seg.left = seg.right - remaining;
			    }
			    intersect = fullworklist.intersection( seg );
			}

			if (trueremaining > 0) {
			    remaining -= intersect;
			    newsegs.add( seg );
			    // now update the remainder counts for the isochores.  Here, overshoot is very
			    // possible.  If that happens, charge a random other isochore.
			    int rounds = 0;
			    while (intersect > 0) {
				long isoleft = (Long)isoremaining.get(iso);
				if (isoleft >= intersect) {
				    isoleft -= intersect;
				    intersect = 0;
				} else {
				    intersect -= isoleft;
				    isoleft = 0;
				}
				isoremaining.set(iso, isoleft);
				iso = (iso + 1) % isoremaining.size();
				rounds += 1;
				if (rounds > isoremaining.size()) {
				    // iso counts are reset later
				    intersect = 0;
				}
			    }
			}
		    }

		    if (dumpSegments) {
			intersectedSegments.dump(segmentKeyString + "\t" + chromosome);
		    }


		    // Finally, simply compute the intersections with the relevant workspaces for each isochore
		    for (idx = 0; idx < isochores.size(); idx++) {
			String isochore = (String)isochores.get(idx);
			currentSegments.put(isochore, null);                       // reduce memory usage
			intersectedSegments = new SegmentList( (SegmentList)workspace.get( isochore ) );
			intersectedSegments.intersect( segList );
			currentWorkSegs.put(isochore, intersectedSegments);

		    }
		}
    }


	int removeLowSupportAnnotations( double densityThreshold ) {

		Set outside = new HashSet();
		Iterator iter = annotationStats.keySet().iterator();
		while (iter.hasNext()) {
			String ann = (String)iter.next();
			Statistics stats = (Statistics)annotationStats.get( ann );
			double avg = stats.average();
			// criterion
			if ((avg<densityThreshold)||(stats.getObserved()<densityThreshold)) {
				outside.add( ann );
			}
		}

		Iterator annIter = outside.iterator();
		while (annIter.hasNext()) {
			String ann = (String)annIter.next();
			annotationStats.remove( ann );
			if (verbose > 2) {
				System.out.println("# Removing low density annotation "+ann);
			}
		}
		return outside.size();

	}


	private double histogramTotal( int[] hist ) {

		double sum = 0.0;
		for (int i=0; i<hist.length; i++) {
			sum += hist[i] * i;
		}
		return sum;
	}	



	private double histogramMean( int[] hist ) {

		double n = 0.0;
		for (int i=0; i<hist.length; i++) {
			n += hist[i];
		}	
		if (n < 1.0) {
			n = 1.0;
		}
		return histogramTotal( hist )/n;
	}

	private double histogramVar( int[] hist ) {

		double n = 0.0;
		double var = 0.0;
		double mean = histogramMean( hist );
		for (int i=0; i<hist.length; i++) {
			var += hist[i] * (i - mean)*(i - mean);
			n += hist[i];
		}
		if (n < 2.0) {
			n = 2.0;
		}
		return var / (n-1.0);
	}



	void dumpResults( double zThreshold, double pThreshold, int numSamples, boolean dumpSamples) {

		FalseDiscoveryStats fds = new FalseDiscoveryStats();

		try {

			if (0.01 <= fdrmaxp) {
				fds.dumpFDRSummary( annotationStats, 0.01, numSamples, zThreshold );
			}	
			if (0.005 <= fdrmaxp) {
				fds.dumpFDRSummary( annotationStats, 0.005, numSamples, zThreshold );
			}
			if (0.002 <= fdrmaxp) {
				fds.dumpFDRSummary( annotationStats, 0.002, numSamples, zThreshold );
			}
			if (0.001 <= fdrmaxp) {
				fds.dumpFDRSummary( annotationStats, 0.001, numSamples, zThreshold );
			}
			if (0.0005 <= fdrmaxp) {
				fds.dumpFDRSummary( annotationStats, 0.0005, numSamples, zThreshold );
			}
			if (0.0002 <= fdrmaxp) {
				fds.dumpFDRSummary( annotationStats, 0.0002, numSamples, zThreshold );
			}
			if (0.0001 <= fdrmaxp) {
				fds.dumpFDRSummary( annotationStats, 0.0001, numSamples, zThreshold );
			}
			if (0.00005 <= fdrmaxp) {
				fds.dumpFDRSummary( annotationStats, 0.00005, numSamples, zThreshold );
			}
			if (0.00002 <= fdrmaxp) {
				fds.dumpFDRSummary( annotationStats, 0.00002, numSamples, zThreshold );
			}
			if (0.00001 <= fdrmaxp) {
				fds.dumpFDRSummary( annotationStats, 0.00001, numSamples, zThreshold );
			}

		} catch (Exception e) {
			System.out.println("# Internal error: "+e.getMessage()+"\n# (Proceeding)");
		}

		DecimalFormat df = new DecimalFormat();

		// Sort annotations according to relative representation
		List annotationList = new ArrayList(0);
		annotationList.addAll( annotationStats.values() );
		Collections.sort( annotationList, new relRepComparator() );

		// Loop over all annotations
		// Loop over all annotations
		System.out.println("Z\tOver/underrepr(%)\tPvalue\tObserved\tExpected\tLower 95CI\tUpper 95CI\tStdDev\tAnnotation");
		Iterator iter = annotationList.iterator();
		while (iter.hasNext()) {
			Statistics stats = (Statistics)iter.next();

			double expPval = stats.experimentalP();
			double relrep = stats.relativeRepresentation();
			double z = stats.z();

			// criterion
			if ((Math.abs(z)>=zThreshold) && (Math.abs(expPval)<=pThreshold)) {

				df.applyPattern("###.##");
				String devexpStr = df.format(z);
				String relrepStr = df.format(relrep);
				df.applyPattern("#.0000000");
				String experpStr = df.format(Math.abs(expPval));
				String propStr = df.format(stats.getObserved());
				String epropStr = df.format(stats.average());
				String ci0 = df.format(stats.low95CI());
				String ci1 = df.format(stats.high95CI());
				String sdStr = df.format(stats.sd());
				System.out.println(devexpStr+"\t"+relrepStr+"\t"+experpStr+"\t"+propStr+"\t"+epropStr+"\t"+ci0+"\t"+ci1+"\t"+sdStr+"\t"+stats.getAnnotation());

			}
		}

		System.out.println("\n");

	}





	class relRepComparator implements Comparator {

		public int compare(Object o1, Object o2) {

			Statistics s1 = (Statistics)o1;
			Statistics s2 = (Statistics)o2;
			double r1 = s1.relativeRepresentation();
			double r2 = s2.relativeRepresentation();
			if (r1-r2 < 0.0) {
				return -1;
			} else if (r1-r2 > 0.0) {
				return 1;
			} 
			return 0;
		}

		public boolean equals(Object o1, Object o2) {

			return compare(o1,o2) == 0;

		}

	}


	void dumpSummary( int originalNumIDs, int numNotInWorkspace, int numPurged, int numEmptyAnnotations, int numDuplicatedAnnotations, int numUndefIDs ) {

		long totalWork = 0;
		long totalSegs = 0;
		long totalSegsWorksp = 0;
		long totalIds = 0;
		long totalIdLength = 0;
		long totalIdWorksp = 0;

		DecimalFormat df = new DecimalFormat();
		df.applyPattern("#,###");

		// Sum segments over base chromosomes, not isochores
		// Note: 'segments' contains both base and isochore keys, mapping to same segment list
		Iterator iter = synonyms.entrySet().iterator();
		boolean bailOut = false;
		while (iter.hasNext()) {

			Entry e = (Entry)iter.next();
			String chrom = (String)e.getKey();
			SegmentList sl = (SegmentList)segments.get( chrom );
			if (sl == null) {
			    System.out.println("No segments found for base chromosome "+chrom+" listed in synonym file -- is it defined in the workspace?");
			    bailOut = true;
			} else {
			    totalSegs += sl.length();
			}

		}

		if (bailOut) {
		    throw new Error("Bailing out because of workspace/synonym mismatch.");
		}

		// Loop over isochores
		if (this.verbose > 2) {
			System.out.println("# Isochore and segment statistics:");
		}
		iter = (new TreeSet(workspace.keySet())).iterator();
		while (iter.hasNext()) {

			String isoch = (String)iter.next();
			SegmentList sl = (SegmentList)workspace.get( isoch );
			totalWork += sl.length();
			long segWorksp = ((Long)totalLength.get(isoch)).longValue();
			if (this.verbose > 2) {
				int[] hist = (int[])lengthDistr.get( isoch );
				double density;
				if (sl.length() >= 1) {
					density = segWorksp / (double)sl.length();
				} else {
					density = 0.0;
				}
				double mean = histogramMean(hist);
				double sd = Math.sqrt( histogramVar(hist) );
				System.out.println("#  "+isoch+":\tsize="+sl.length()+" segment density="+density+" mean length="+mean+" sd="+sd);
			}
			totalSegsWorksp += segWorksp;
		}

		iter = ids.entrySet().iterator();
		HashMap id2worksp = new HashMap(0);
		while (iter.hasNext()) {

			totalIds += 1;
			Entry e = (Entry)iter.next();
			Location loc = ((Location)(e.getValue()));
			// Actual size of segment
			totalIdLength += loc.getLength();
			// Compute size of intersection with workspace, on each isochore
			long intersection = 0;
			List syns = (List)synonyms.get( loc.chrom ); 
			for (int idx=0; idx<syns.size(); idx++) {
				String isoch = (String)syns.get(idx);
				SegmentList segList = (SegmentList)workspace.get( isoch );
				if (segList != null) {
					intersection += segList.intersection( loc.getSegments() );
				}
			}

			totalIdWorksp += intersection;
			id2worksp.put( e.getKey(), new Long(intersection) );
		}

		// Calculate statistics for annotation
		iter = annotation.entrySet().iterator();
		List annworksp = new ArrayList(0);
		long grandtotal = 0;
		while (iter.hasNext()) {
			Entry e = (Entry)iter.next();
			long[] ids = (long[])e.getValue();
			long total = 0;
			for (int i=0; i<ids.length; i++) {
				Long id = new Long( ids[i] );
				Long nucs;
				if (id2worksp.containsKey(id)) {
					nucs = (Long)id2worksp.get( id );
				} else {
					nucs = new Long(0);
				}
				total += nucs.longValue();
			}
			grandtotal += total;
			annworksp.add( new Long(total) );
		}
		if (annworksp.size() == 0) {
			throw new Error("Error - no annotations defined.");
		}
		Collections.sort( annworksp );
		long median = ((Long)annworksp.get( annworksp.size()/2 )).longValue();
		long tenperc = ((Long)annworksp.get( annworksp.size()/10 )).longValue();
		long ninetyperc = ((Long)annworksp.get( (9*annworksp.size())/10 )).longValue();
		long mean = grandtotal / annworksp.size();

		System.out.println("# Removed "+numPurged+" IDs not referred to in annotations");
		System.out.println("# Removed "+numNotInWorkspace+" IDs not overlapping workspace");
		System.out.println("# Found "+numUndefIDs+" unique undefined IDs referred to in annotations");
		System.out.println("# Removed "+numEmptyAnnotations+" empty annotations");
		System.out.println("# Removed "+numDuplicatedAnnotations+" duplicated annotations");
		System.out.println("# Workspace:\t"+df.format(totalWork)+" nucs");
		System.out.println("# Ids:      \t"+totalIds+" (before purge: "+originalNumIDs+") covering "+df.format(totalIdLength)+" nucs of which "+df.format(totalIdWorksp)+" in workspace");
		System.out.println("# Segments:\t"+df.format(totalSegs)+" nucs of which "+df.format(totalSegsWorksp)+" in workspace");
		System.out.println("# Annotations:\t"+annotation.size()+"; mean overlap "+df.format(mean)+", median "+df.format(median)+", 10%-90% quantiles "+df.format(tenperc)+"-"+df.format(ninetyperc));
	}


	// Input data

    Map workspace;       // String (chromosome/isochore) -> SegmentList
    Map fullWorkspaces;   // String (chromosome) -> SegmentList; used for unintersected sampling
	Map segments;        // String (chromosome) -> SegmentList
	Map ids;             // Long (id) -> Location
	Map annotation;      // String (annotation) -> array of longs
	Map original;        // Maps chromosome synonym to base chromosome
	Map synonyms;        // Maps base chromosome to list of synonyms
	Map subsets;         // Maps subset id to list of synonyms (default: "" to all)
	Random random;

	// Derived data

	//long workIgsLength;  // Total length of IGS intersecting workspace
	Map workIgsLength;   // Total length of IGS intersecting workspace within subset
	Map sampleIgsLength; // Total length of externally IGS intersecting workspace within subset;
	Map totalLength;     // String (chromosome) -> Long (total segments length intersecting with workspace)
	Map lengthDistr;     // String (chromosome) -> int[] (length distribution of segments intersecting w. worksp.)
	Map largeWorkspace;  // String (chromosome) -> SegmentList; enlarged workspace
	Map currentSegments; // Current set of (sampled) segments; == segments at first iteration
	Map currentWorkSegs; // currentSegments intersected with workspace
	Map idIntersection;  // Long (id) -> Integer, intersection with workspace and currentSegments
	Map annotationStats; // String (annotation) -> Statistics
	Map segmentListSamplers;   // Initialization is expensive and depends only on largeWorkspace


}

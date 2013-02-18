package segment;


import java.util.List;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Iterator;
import java.text.ParsePosition;

import app.Parser;
import app.Globals;
import segment.Segment;

/***
 * 
 * List of non-intersecting segments
 * 
 */


public class SegmentList {
	
	/**
	 * 
	 * Creates a segment list from a string consisting of tab-delimited records of the form
	 * (start,end).  Whitespace is skipped.
	 * 
	 * @param input  string of form "(start,end)\t(start,end)\t...."
	 * @param offset offset in string to begin parsing
	 */
	public void parseSegmentList( String input, int offset, int shoulder ) {

		// Create the list
		segList = new ArrayList(0);
		
		// Parses the string from 'offset'.
		// The string is expected to consist of segments of the form "(start,end)",
		//  with \t separating segments.  Whitespace is skipped
		
		ParsePosition parsepos = new ParsePosition(0);

		parsepos.setIndex( offset );
		while (parsepos.getIndex() < input.length()) {

			try {
				
				Segment s = Segment.parse( input, parsepos );
				if (s==null) {
					throw new Error("Out of memory!");
				}
				if (shoulder == 0) {
				    segList.add( s );
				} else {
				    if (s.length() > 0) {
					segList.add( new Segment( s.left()-shoulder, s.right()+shoulder ) );
				    }
				}
				
			} catch (Segment.SegmentParseDone e) {
				break;
			}
			
		}
		
		// Check that we're at end-of-line
		Parser.skipWhitespace( input, parsepos );
		if (input.length() < parsepos.getIndex() ) {
		    throw new Error("Parse problem - line does not end in newline at position "+parsepos.getIndex());
		}
		
		// Sort, merge, and remove empties
		normalize();
	}


	public SegmentList( String input, int offset, int shoulder ) {

	    parseSegmentList( input, offset, shoulder );

	}


	public SegmentList( String input, int offset ) {

	    parseSegmentList( input, offset, 0 );

	}


	public SegmentList( SegmentList segments ) {

		segList = new ArrayList( segments.segList.size() );
		segList.addAll( segments.segList );
		
		checkConsistency();
		
	}
	
	
	public SegmentList() {
		
		segList = new ArrayList(0);
		
	}


        /* Single-segment segmentlist */
        public SegmentList( Segment seg ) {
	    
	    segList = new ArrayList(1);
	    segList.add( seg );
	    
	}

	
	private void normalize() {
	
		// get array underlying the list
		Segment[] array = (Segment[])segList.toArray( new Segment[0] );
		
		// sort it in ascending order of 'left' field
		Arrays.sort( array );
		
		// merge
		int i=0;
		while (i < array.length && array[i].length()==0) {
			i++;
		}
		if (i < array.length) {
			Segment curseg = array[i];
			int cursegidx = 0;
			for (; i<array.length; i++) {
				Segment newseg = array[i];
				if (newseg.length() > 0) {
					if (newseg.left() <= curseg.right() ) {
						// Overlap - extend current segment to include new one
						curseg = new Segment( curseg.left(), Math.max( curseg.right(), newseg.right() ) );
						// Remove new segment
						array[i] = Segment.emptySegment();
					} else {
						// No overlap - store current segment, and set new to current
						array[cursegidx] = curseg;
						curseg = newseg;
						cursegidx = i;
					}
				}
			}
			// Store final segment
			array[cursegidx] = curseg;
		}
		
		// Remove empty segments
		int cursegidx = 0;
		for (i=0; i<array.length; i++) {
			if (array[i].length() != 0) {
				array[cursegidx] = array[i];
				cursegidx += 1;
			}
		}
		
		// Put into segList again
		segList.clear();
		segList.addAll( (Arrays.asList( array )).subList(0,cursegidx) );
		
		checkConsistency();

	}
	
	
	/**
	 * 
	 * @param extend  size by which each segment is extended (in both directions)
	 */
	public void extend( long sizeleft, long sizeright ) {
		
		// We rely on the contract that the segmentList does not contain empty segments
		for (int i=0; i<segList.size(); i++) {
			
			Segment s = (Segment)segList.get(i);
			segList.set(i, new Segment( s.left() - sizeleft, s.right() + sizeright ) );
			
		}
		
		// Merge
		normalize();
		
	}


	/**
	 * Returns list of segments that intersect (partially or completely) with those in other.

	 * @param other SegmentList
	 * @return indices of this.segList that intersect other
	 */
	public List intersectingSegments( SegmentList other ) {

		int selfidx = 0;
		int otheridx = 0;
		List segments = new ArrayList(0);
	
		while (selfidx < segList.size() && otheridx < other.segList.size() ) {
		
			// Deal with all other segments with right endpoint <= that of self segment
			boolean intersects = false;
			while (otheridx < other.segList.size() && ((Segment)other.segList.get(otheridx)).right() <= ((Segment)segList.get(selfidx)).right() ) {

				// Calculate intersection
				Segment intersection = ((Segment)segList.get(selfidx)).intersection( ((Segment)other.segList.get(otheridx)));
				if (intersection.length() > 0) {
					intersects = true;
				}
			
				otheridx += 1;
			
			}

			// Get intersection of self segment with other,
			if (otheridx < other.segList.size()) {
				Segment intersection = ((Segment)segList.get(selfidx)).intersection( ((Segment)other.segList.get(otheridx)));
				if (intersection.length() > 0) {
					intersects = true;
				}
			}
			
			if (intersects) {
				segments.add( segList.get(selfidx) );				
			}
			
			selfidx += 1;
		}

		return segments;
	}

	
	
	/**
	 * 
	 * Intersect this segment list with another
	 * 
	 * @param segmentlist
	 */
	public void intersect( SegmentList other ) {
		
		ArrayList newList = new ArrayList(0);
		
		int selfidx = 0;
		int otheridx = 0;
		int selfsize = segList.size();
		int othersize = other.segList.size();
		Segment nullseg = Segment.emptySegment();
		List othersegList = other.segList;                     // for speed
		
		Segment otherSeg = Segment.emptySegment();             // to make the compiler happy
		while (selfidx < selfsize && otheridx < othersize ) {
			
			Segment selfSeg = (Segment)segList.get(selfidx);
			if (otheridx < othersize) {
				otherSeg = (Segment)othersegList.get(otheridx);
				// Deal with all other segments with right endpoint <= that of self segment
				while (otherSeg.right <= selfSeg.right ) {

					// Calculate intersection; use that the empty segment == nullseg
					Segment intersection = selfSeg.intersection( otherSeg );
					if (intersection != nullseg ) {
						newList.add( intersection );
					}
				
					otheridx += 1;
					if (otheridx < othersize) {
						otherSeg = (Segment)othersegList.get(otheridx);
					} else {
						break;
					}
				}
			}
			
			// Get intersection of self segment with other,
			if (otheridx < othersize) {
				Segment intersection = selfSeg.intersection( otherSeg );
				if (intersection != nullseg) {
					newList.add( intersection );
				}
			}
			
			selfidx += 1;
		}
		
		// normalization is not necessary
		segList = newList;

		checkConsistency();
		
	}
	
	
	/**
	 * 
	 * Take the union of this segment list with another
	 * 
	 * @param segmentlist
	 */
	public void union( SegmentList segmentlist ) {
		
		addList( segmentlist.segList );
		
	}
	
	
	public void addList( List segs ) {

		segList.addAll( segs );
		normalize();
		checkConsistency();		

	}
	

	/**
	 * Compute the intersection of segment list with a single segment
	 * 
	 * @param seg
	 * @return total number of sites overlapping with seg
	 */
	public long intersection( Segment seg ) {

		int idx = Collections.binarySearch( segList, seg );
		
		if (idx<0) {
			idx = -(idx+1);       // insertion point
			if (idx>0) idx -= 1;  // left neighbour could overlap seg
		}
		
		int sls = segList.size();
		if (idx>=sls) {
			return 0;
		}
		
		long count = 0;
		Segment curseg = (Segment)segList.get(idx);
		do {
			count += ((curseg).intersection( seg )).length();
			idx += 1;
			if (idx<sls) {
				curseg = (Segment)segList.get(idx);
			} else {
				break;
			}
		} while (curseg.left <= seg.right);
		
		return count;
		
	}


    /** Trims segment list, by removing nucleotides starting from the segment that includes pos */

    public void trim( long pos, long size ) {

	Segment seg = new Segment(pos,pos+1);
	int idx = Collections.binarySearch( segList, seg );
	if (idx<0) {
	    idx = -(idx+1);       // insertion point
	    if (idx == segList.size()) idx=0;
	}
	while (size > 0) {
	    seg = (Segment)segList.get(idx);
	    if (seg.length() < size) {
		segList.set(idx, new Segment(0,0) );
		size -= seg.length();
	    } else {
		segList.set(idx, new Segment( seg.left() + size, seg.right() ) );
		size = 0;
	    }
	    idx += 1;
	    if (idx == segList.size())
		idx = 0;
	}
    }



        /** 
	 * Compute the length of the intersection with another SegmentList.  This implementation is
	 * sluggish for long dense arrays, because of the repeated binary searches; the ordering of
	 * the SegmentList could be exploited to speed things up here.
	 */
        public long intersection( SegmentList segs ) {

	    Iterator iter = segs.segList.iterator();
	    long length = 0;

	    while (iter.hasNext()) {

		Segment seg = (Segment)iter.next();
		length += intersection( seg );

	    }

	    return length;

	}
		
	
	public long length() {
		
		return intersection( new Segment( Integer.MIN_VALUE, Integer.MAX_VALUE ) );
		
	}
	
	
	public long max() {
		
		return ((Segment)segList.get( segList.size() - 1 )).right();

	}
	
	
	public int numSegments() {
		return segList.size();
	}
	
	
	public Segment get(int i) {
		return (Segment)segList.get(i);
	}


    public void dump( String s ) {

	System.out.print(s);
	for (int i=0; i<segList.size(); i++) {
	    
	    System.out.print("\t");
	    ((Segment)segList.get(i)).dump();

	}
	System.out.println("");
    }
	    
	    
	

	public void checkConsistency() {
		
		if (!Globals.doChecks) {
			return;
		}
		// Check for: no empties, no overlaps (but touching is allowed), and no HUGE segments
		for (int i=0; i<segList.size(); i++) {
			
			Segment s = (Segment)segList.get(i);
			if (s.left() >= s.right()) {
				throw new Error("Inconsistent segList - empty segment");
			}
			if (s.length() > 1000000000000L) {
				throw new Error("Inconsistent segList - segment over 1Tb");
			}
			if (i>0) {
				Segment sl = (Segment)segList.get(i-1);
				if (sl.right() == s.left()) {
					throw new Error("Inconsistent segList - segments touch");
				}
				if (sl.right() > s.left()) {
					throw new Error("Inconsistent segList - segments overlap");
				}
			}
			
		}
		
	}
	
	private List segList;
	
}


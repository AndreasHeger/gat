package app;

import segment.Segment;
import segment.SegmentList;



public class Location {
	
    Location( String chrom, Segment seg ) {
	
	this.chrom = chrom.intern();
	this.seg = seg;
	
    }
    
    Location( String chrom, SegmentList segments ) {
	
	this.chrom = chrom.intern();
	this.seg = segments;
	
    }
    
    
    SegmentList getSegments() {
	
	if (seg instanceof SegmentList) {
	    
	    return (SegmentList)seg;
	    
	} else {

	    /* Create single-element SegmentList */	    
	    SegmentList list = new SegmentList( (Segment)seg );
	    return list;
	    
	}
	
    }
    
    
    long getLength() {
	
	if (seg instanceof Segment)
	    return ((Segment)seg).length();

	return ((SegmentList)seg).length();

    }


    String chrom;
    Object seg;     /* a Segment or a SegmentList */

    
}

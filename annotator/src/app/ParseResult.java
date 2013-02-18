package app;
	
import segment.Segment;
import segment.SegmentList;
import java.util.List;

public class ParseResult {
		
	final static int NONE = -1;
	final static int WORK = 0;
	final static int SEGS = 1;
	final static int ANN = 2;
	final static int ID = 3;
	final static int SYN = 4;
	final static int SUBSET = 5;
	
	/**
	 * Constructor for WORK and SEGS
	 */
	ParseResult( int type, String chr, SegmentList sl ) {
		
		this.type = type;
		this.chr = chr;
		this.seglist = sl;
		
	}

	/**
	 * Constructor for ANN
	 */
	ParseResult( int type, String annotation, List list ) {
		
		this.type = type;
		this.ann = annotation;
		this.list = list;
		
	}

	/**
	 * Constructor for ID
	 */
	ParseResult( long id, String chromosome, Segment seg ) {
		
		this.type = ID;
		this.chr = chromosome;
		this.id = id;
		this.seg = seg;
		
	}

	/**
	 * Constructor for SYN and SUBSET
	 *
	 */
	ParseResult( int type, String synonym, String chromosome ) {

		this.type = type;
		this.chr = chromosome;    // base (syn) or subset id (subset)
		this.ann = synonym;       // synonym
	}
	
	
	ParseResult() {
		
		type = NONE;
		
	}


    void setIDSegList( long id, String chromosome, SegmentList list ) {

	    this.type = ID;
	    this.chr = chromosome;
	    this.id = id;
	    this.seg = null;
	    this.seglist = list;
	
	}
	
	int type;
	String chr;
    String ann;
	List list;
	SegmentList seglist;
	long id;
	Segment seg;
	
}

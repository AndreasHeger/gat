package segment;


import java.text.DecimalFormat;
import java.text.ParsePosition;

import app.Parser;


public final class Segment implements Comparable {

	public Segment( long left, long right ) {
		
		this.left = left;
		this.right = right;
		
	}
	
	public final long length() {
		
		if (left<right) {
			return right-left;
		} else {
			return 0;
		}
	}
	
	public final long left() {
		return this.left;
	}
	
	public final long right() {
		return this.right;
	}
	
	public final int compareTo( Object o ) {
	
		Segment s = (Segment)o;
		if (left < s.left) {
			return -1;
		} else if (left > s.left ) {
			return 1;
		}
		return 0;
	
	}

        public final void dump () {
	    System.out.print("("+this.left+","+this.right+")");
	}
	

	
	public final Segment intersection( Segment s ) {
		
		long l = s.left;
		long r = s.right;
		if (left > l) l=left;
		if (right < r) r=right;
		if (l >= r) return nullSegment;
		if (l==left && r==right) return this;
		if (l==s.left && r==s.right) return s;
		return new Segment( l, r );
			
	}
	
	public final Segment union( Segment s ) {
		
		if (!touch(s)) {
			throw new Error("Un-unionable");
		}
		return new Segment( Math.min(left,s.left), Math.max(s.right, right) );
	}

	
	public final boolean touch( Segment s ) {
		
		return Math.max(left,s.left) <= Math.min(right,s.right);		
		
	}
	
	static final public Segment emptySegment() {
		
		return nullSegment;
		
	}

	static public final Segment parse( String input, ParsePosition parsepos ) throws SegmentParseDone {

		DecimalFormat numberparser = new DecimalFormat();
		numberparser.setGroupingUsed( false );

		int startpos = parsepos.getIndex();
		Parser.skipWhitespace(input,parsepos);
		if (input.charAt(parsepos.getIndex()) != '(') {
			throw new SegmentParseDone();
		}
		Parser.skipWhitespace( input, parsepos, 1 );
		Number start = numberparser.parse( input, parsepos );
		if (start == null) {
			throw new Error("Parse problem at position "+parsepos.getIndex());
		}
		Parser.skipWhitespace(input, parsepos);
		if (input.charAt(parsepos.getIndex()) == 'L') {
		    Parser.skipWhitespace(input, parsepos, 1);
		}
		if (input.charAt(parsepos.getIndex()) != ',') {
			throw new Error("Parse problem (no comma) at position "+parsepos.getIndex());
		}
		Parser.skipWhitespace( input, parsepos, 1 );
		Number end = numberparser.parse( input, parsepos );
		if (end == null) {
			throw new Error("Parse problem (second number) at position "+parsepos.getIndex());
		}
		Parser.skipWhitespace( input, parsepos );
		if (input.charAt(parsepos.getIndex()) == 'L') {
		    Parser.skipWhitespace(input, parsepos, 1);
		}
		if (input.charAt(parsepos.getIndex()) != ')') {
			throw new Error("Parse problem (no closing bracket) at position "+parsepos.getIndex());
		}
		Parser.skipWhitespace(input, parsepos, 1 );
		
		long s = start.longValue();
		long e = end.longValue();
		if (s>=e) {
		    //			throw new Error("Segment '"+input.substring(startpos, parsepos.getIndex())+"' parsed to empty seg s="+s+", e="+e+" at position "+parsepos.getIndex());
		    System.out.println(" Warning: Segment '"+input.substring(startpos, parsepos.getIndex())+"' parsed to empty seg s="+s+", e="+e+" at position "+parsepos.getIndex());
		}
		long maxsize = 1000000000000L;  // 1 Tb
		if ((e-s > maxsize) || (e > maxsize) || (s < -maxsize)) {
			throw new Error("Segment '"+input.substring(startpos, parsepos.getIndex())+"' parsed to oversize seg s="+s+", e="+e+" at position "+parsepos.getIndex());			
		}
		
		return new Segment( start.longValue(), end.longValue() );
		
	}

	public static final class SegmentParseDone extends Exception {
		static final long serialVersionUID = 0;
	}
	
	public long left;
	public long right;
	
	private static final Segment nullSegment = new Segment(0,0);
	
}

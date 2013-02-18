package app;


import java.text.ParsePosition;
import java.util.List;
import java.util.ArrayList;
import java.io.LineNumberReader;
import java.io.IOException;

import app.ParseResult;
import segment.Segment;
import segment.SegmentList;

public class Parser {
	
	
	Parser(AnnotationData ad) {
		
		this.ad = ad;
		worklines = annlines = synlines = seglines = idlines = 0;
		
	}
	
	
	void parseAnnotationData( LineNumberReader reader, int verbose, int shoulder ) {
		
		try {
			
			String line = reader.readLine();
			while (line != null) {
				
				ParseResult pr;
				try {
					pr = parse(line,reader.getLineNumber(),shoulder);
				} catch (Error e) {
					throw new Error(e.getMessage()+" at line "+reader.getLineNumber());
				}
				if (pr == null) {
					throw new Error("Out of memory!  Use the -Xmx option to increase the maximum Java heap size.");
				}
				switch (pr.type) {
				case ParseResult.NONE:
					break;
				case ParseResult.WORK:
				    if (ad.workspace.containsKey( pr.chr )) {
					pr.seglist.intersect( (SegmentList)ad.workspace.get( pr.chr ) );
				    }
					ad.workspace.put( pr.chr, pr.seglist );
					break;
				case ParseResult.SEGS:
					ad.segments.put( pr.chr, pr.seglist );
					break;
				case ParseResult.ANN:
				    long[] annlist = new long[ pr.list.size() ];
				    for (int i=0; i<pr.list.size(); i++) {
					annlist[i] = ((Long)pr.list.get(i)).longValue();
				    }
				    ad.annotation.put( pr.ann, annlist );
				    break;
				case ParseResult.ID:
   				        if (pr.seg != null) {
					    ad.putIds( new Long(pr.id), new Location( pr.chr, pr.seg ) );
					} else {
					    ad.putIds( new Long(pr.id), new Location( pr.chr, pr.seglist ) );
					}
					break;
				case ParseResult.SYN:
					ad.original.put( pr.ann, pr.chr );    // pr.ann = derived, pr.chr = original
					if ( !ad.synonyms.containsKey( pr.chr )) {
					    ad.synonyms.put( pr.chr, new ArrayList(0) );
					}
					if ( !((List)ad.synonyms.get( pr.chr )).contains( pr.ann ) ) {
						((List)ad.synonyms.get( pr.chr )).add( pr.ann );
					}
					break;
				case ParseResult.SUBSET:
					if ( !ad.subsets.containsKey( pr.chr )) {
						ad.subsets.put( pr.chr, new FastArrayList(0) );
					}
					if ( !((List)ad.subsets.get( pr.chr )).contains( pr.ann ) ) {
						((List)ad.subsets.get( pr.chr )).add( pr.ann );
					}
					break;
				default:
					throw new Error("Unknown parse result obtained ("+pr.type+")");
				}

				line = reader.readLine();
				
			}
			
		} catch (IOException io) {
			
			throw new Error(io.getMessage());
			
		}
	}
	


	public final ParseResult parse( String line, int lineno, int shoulder ) {
		
		if (line.length() < 3) {
			return new ParseResult();   // none
		}
		
		if (line.charAt(0) != '#' || line.charAt(1) != '#') {
			return new ParseResult();   // none
		}
				
		int type = -1;
		ParsePosition pp = new ParsePosition( 2 );
		String keyword = parseWord( line, pp );
		
		if (keyword.equals("Work")) {
			type = ParseResult.WORK;
			worklines += 1;
		} else if (keyword.equals("Segs") || keyword.equals("Seg")) {
			type = ParseResult.SEGS;
			seglines += 1;
		} else if (keyword.equals("Ann")) {
			type = ParseResult.ANN;
			annlines += 1;
		} else if (keyword.equals("Id")) {
			type = ParseResult.ID;
			idlines += 1;
		} else if (keyword.equals("Synonym")) {
			type = ParseResult.SYN;
			synlines += 1;
		} else if (keyword.equals("Subset")) {
			type = ParseResult.SUBSET;
			subsetlines += 1;
		} else {
			throw new Error("Error parsing - keyword '"+keyword+"' not recognized");
		}
		
		
		// First parameter is either chromosome or annotation or subset id.  Intern the string for memory efficiency
		String first = parseWord( line, pp ).intern();

		// Work and Segs: chr, segs
		if (type == ParseResult.WORK || type == ParseResult.SEGS) {
			
		    SegmentList sl = null;
		    if (type == ParseResult.SEGS) {
		        sl = new SegmentList( line, pp.getIndex(), shoulder );
		    } else {
			sl = new SegmentList( line, pp.getIndex() );
		    }
		    return new ParseResult( type, first, sl );
			
		}
		
		// Synonym or subset: chr, chr
		if (type == ParseResult.SYN || type == ParseResult.SUBSET) {
			
			String chrom = parseWord( line, pp, true ).intern();
			return new ParseResult( type, first, chrom );
			
		}
		
		// Ann: annotation, ids
		if (type == ParseResult.ANN) {
		
			List idl = new ArrayList(0);
			String word = parseWord( line, pp, true );
			while (word.length() > 0) {
			
				try {
					
					long id = Long.parseLong( word );
					idl.add( new Long(id) );
					word = parseWord( line, pp, true );
					
				} catch (NumberFormatException e) {
					
					word = "";
					System.out.println("Possible problem at position "+pp.getIndex()+" in line "+lineno+" defining annotation "+first);
					
				}
			}
			return new ParseResult( type, first, idl );
		}

		// Id: id, chr, seg(s?)
		String chrom = parseWord( line, pp );
		Segment seg = null;
		SegmentList sl = null;
		try {
		    sl = new SegmentList( line, pp.getIndex() );
		    if (sl.numSegments() == 1) {
			// If just one, store as a single segment
		    	seg = sl.get(0);
		    }
		    if (sl.numSegments() == 0) {
		    	throw new Error("No segments...");
		    }
		} catch (Error e) {
			// Allow problems here - ID file is a bit dodgy
			System.out.println("Warning -- segment parse error at line "+lineno+" (ID="+first+" chrom="+chrom+"): "+e.getMessage());
			seg = new Segment(0,0);
		}
		if (seg != null) {
		    // Return a single-segment ID
		    return new ParseResult( Long.parseLong(first), chrom, seg);
		} else {
		    // Return a segment-list ID
		    ParseResult pr = new ParseResult();
		    pr.setIDSegList( Long.parseLong(first), chrom, sl );
		    return pr;
		}

	}
	
	
	public static final void skipWhitespace( String input,  ParsePosition parsepos, int skip ) {
		// skip whitespace
		int pointer = parsepos.getIndex() + skip;
		while (pointer < input.length() && (input.charAt(pointer) == ' ' || input.charAt(pointer) == '\t')) {
			pointer += 1;
		}
		parsepos.setIndex( pointer );
	}

	
	public static final void skipWhitespace( String input, ParsePosition parsepos ) {
		Parser.skipWhitespace( input, parsepos, 0 );
	}
	
	
	public static final String parseWord( String input, ParsePosition pp, boolean delimitBySpace ) {
//hotspot 9%		
		// Skip any leading spaces
		skipWhitespace( input, pp );
		// Word extends until next \t (or space)
		int p = pp.getIndex();
		while (p < input.length() && input.charAt(p) != '\t' && !(delimitBySpace && input.charAt(p) == ' ')) {
			p++;
		}
		String result = input.substring(pp.getIndex(), p);
		result = result.trim();
		pp.setIndex( p );
		// Skip \t (and any whitespace of next word)
		skipWhitespace( input, pp );
		
		return result;
		
	}
	
	public static final String parseWord( String input, ParsePosition pp) {
		return parseWord( input, pp, false );
	}

    public void dumpSummary(int verbose) {
	
    	if (verbose>0) {
    	
    		System.out.println("# Parsed "+seglines+" segment lines");
    		System.out.println("# Parsed "+idlines+" id lines");
    		System.out.println("# Parsed "+annlines+" annotation lines");
    		System.out.println("# Parsed "+worklines+" workspace lines");
    		System.out.println("# Parsed "+synlines+" synonym lines");
    		System.out.println("# Parsed "+subsetlines+" subset lines");
    		
    	}

    	if (seglines == 0) {
    		throw new Error("Parser: no segments definitions found");
    	}
    	if (worklines == 0) {
    		throw new Error("Parser: no workspace definitions found");
    	}
    	if (annlines == 0) {
    		throw new Error("Parser: no annotation lines found");
    	}
    	if (idlines == 0) {
    		throw new Error("Parser: no id definitions found");
    	}
    	
	}

	AnnotationData ad;
        int worklines, seglines, idlines, annlines, synlines, subsetlines;
	
}

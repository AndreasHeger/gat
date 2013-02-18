package app;


import java.util.HashMap;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Iterator;
import java.util.List;
import java.io.LineNumberReader;
import java.io.FileReader;
import java.io.IOException;
import java.io.FileNotFoundException;


import app.AnnotationData;
import segment.SegmentList;


public class Samples {


    String filename;
    AnnotationData ad;
    Parser parser;
    LineNumberReader in;
    HashMap empties;


    Samples(String filename) {

	this.filename = filename;
	this.ad = new AnnotationData();
	this.parser = new Parser( ad );
	this.empties = new HashMap(0);

	try {
	    this.in = new LineNumberReader(new FileReader( filename ));
	} catch (FileNotFoundException e) {	
	    throw new Error("File "+filename+" not found");
	}
    }

    Map getSamples() {
	
	ad.segments = new HashMap();

	String file;
	try {
	    file = in.readLine();
	} catch (IOException e) {
	    throw new Error("Read error reading "+filename);
	}
	if (file == null) {
	    // assume end-of-file
	    return null;
	}

	LineNumberReader segsfile;
	System.out.println("# Reading from file " + file );
	try {
	    segsfile = new LineNumberReader( new FileReader( file ) );
	} catch (FileNotFoundException e) {
	    throw new Error("File "+file+" not found");
	}

	try {
	    parser.parseAnnotationData( segsfile, 1, 0 );
	} catch (Error e) {
	    throw new Error("Error parsing file "+file+":\n"+e.getMessage());
	}

	return ad.segments;

    }

    boolean setSamples( AnnotationData data ) {

	Map samples = getSamples();
	if (samples == null)
	    return false;

	data.sampleIgsLength = new HashMap(0);         // keep lengths here, to be able to compare to originals
	data.currentSegments = new HashMap(0);              

	// Reset total IGS intersection for all subsets
	Iterator subsetiter = data.subsets.keySet().iterator();
	while (subsetiter.hasNext()) {
	    String subset = (String)subsetiter.next();
	    data.sampleIgsLength.put(subset, new Long(0));
	}

	// Loop over all isochores
	Iterator iter = data.workspace.entrySet().iterator();
	while (iter.hasNext()) {
	    Entry e = (Entry)iter.next();
	    SegmentList workList = (SegmentList)e.getValue();
	    String chrom = (String)e.getKey();
	    String basechrom = (String)data.original.get( chrom );

	    SegmentList segList = (SegmentList)samples.get(basechrom);
	    if (segList == null) {
		int num = 1;
		if (empties.get( basechrom ) != null) {
		    num = ((Integer)empties.get(basechrom)).intValue() + 1;
		}
		if (num == 5) {
		    System.out.println("# Warning: No samples for chromosome "+basechrom+" found.  Assuming empty. (Last warning.)");
		} else if (num < 5) {
		    System.out.println("# Warning: No samples for chromosome "+basechrom+" found.  Assuming empty.");
		}
		empties.put( basechrom, new Integer(num) );
		segList = new SegmentList();
		data.currentSegments.put( basechrom, segList );
	    }
	    data.currentSegments.put( chrom, segList );  //!

	    SegmentList currentWorkSeg = new SegmentList( workList );
	    currentWorkSeg.intersect( segList );
	    data.currentWorkSegs.put( chrom, currentWorkSeg );
	
	    // Update sampleIgsLength
	    subsetiter = data.subsets.keySet().iterator();
	    while (subsetiter.hasNext()) {
		String subset = (String)subsetiter.next();
		List syns = (List)data.subsets.get( subset );
		if (syns.contains( chrom )) {
		    // Current isochore 'chrom' is contained within 'subset'
		    long intersection = ((Long)data.sampleIgsLength.get( subset )).longValue() + currentWorkSeg.length();
		    data.sampleIgsLength.put( subset, new Long(intersection) );
		}
	    }
	}

	// Reset total IGS intersection for all subsets
	subsetiter = data.subsets.keySet().iterator();
	while (subsetiter.hasNext()) {
	    String subset = (String)subsetiter.next();
	    long intersection = ((Long)(data.sampleIgsLength.get( subset ))).longValue();
	    long orig = ((Long)(data.workIgsLength.get( subset ))).longValue();
	    double ratio = (intersection + 1.0) / (orig + 1.0);
	    if (ratio < 0.95 || ratio > 1.05) {
		System.out.println("Warning: sampled segment density off: subset="+subset+" original="+orig+" sampled="+intersection);
	    }
	}
	
	return true;

    }
}

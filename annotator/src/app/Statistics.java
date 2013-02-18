package app;


import java.util.Collections;
import java.util.List;
import java.util.ArrayList;
import java.util.Comparator;


/**
 * @author gerton
 *
 */


final class Sample {

    public double prop;    // the sampled proportion
    public int id;         // identifier for this sample

    Sample( double prop, int id ) {
    	this.prop = prop;
    	this.id = id;
    }
}


final class CompareLowSample implements Comparator {
	
	public int compare( Object d, Object e ) {
		double c = ((Sample)d).prop - ((Sample)e).prop;
		if (c<0) {
		    return -1;
		} else if (c>0) {
		    return 1;
		}
		return 0;
	}
}


final class CompareHiSample implements Comparator {
	
	public int compare( Object d, Object e ) {
		double c = ((Sample)d).prop - ((Sample)e).prop;
		if (c<0) {
		    return 1;
		} else if (c>0) {
		    return -1;
		}
		return 0;
	}
}



public final class Statistics {

	private static final double epsilon = 1.0e-15;
	private static final int keepover = 20;

    private static double baseline = 1.0;
    private static int keep = 0;                      // number of low/high samples to keep
    private static int keepLow = 0;                   // same as 'keep' (default), or number of iterations (when dumpSamples is true)
	
	Statistics( String annotation ) {
		
		this.annotation = annotation;
		this.hiSamples = new ArrayList(keep+keepover);
		this.lowSamples = new ArrayList(keepLow+keepover);
		this.highestlow = 1.01;
		this.lowesthigh = -0.01;
		numAbove = 0;
		numBelow = 0;
		num = 0;
		sum = 0.0;
		sumsq = 0.0;
		
	}
	
	static void setBaseline( double baseline ) {
		Statistics.baseline = baseline;
	}
	
	static void setKeep( int keep ) {
		Statistics.keep = keep;
		Statistics.keepLow = keep;
	}
	
	static int getKeep() {
		return Statistics.keep;
	}

    static void keepAll( int numIterations ) {
	Statistics.keepLow = numIterations;
    }

	
	/**
	 * Adds a sample to the set.  The first sample is assumed to be the measurement.
	 * @param proportion
	 */
	final void addSample( double proportion, int id ) {
		
		if (id < 0) {
			observed = proportion;
		} else {
			num += 1;
			sum += proportion;
			sumsq += proportion * proportion;
			if (proportion * baseline >= observed) {
				numAbove += 1;
			}
			// Inclusive comparison is important here; otherwise annotations that
			// overlap with 0 input segments get maximally significant p values...
			if (proportion * baseline <= observed) {
				numBelow += 1;
			}
		}
		// Do not add the observed data (id==-1), and don't add anything if nothing will be kept
		if (id >= 0 && keepLow > 0) {
			Sample s = null;
			if (proportion < highestlow) {
				s = new Sample( proportion, id );
				lowSamples.add(s);
			}
			if (proportion > lowesthigh) {
				if (s == null) {
					s = new Sample( proportion, id );
				}
				hiSamples.add(s);
			}
			if (hiSamples.size() >= keep+keepover || lowSamples.size() >= keepLow+keepover) {
				purge();
			}
		}
	}
	
	
	double average() {
		
		return sum/num;
		
	}
	

	double variance() {
		
		return (sumsq - sum*sum/num) / (num-1);
		
	}
	
	
	double sd() {
		
		return Math.sqrt( variance() );
		
	}
	
	
	double z() {
		
		return (observed - average()) / (sd()+epsilon);
		
	}
	
	
	double low95CI() {
		
		purge();
		if (keepLow >= (int)(0.025 * num) && (int)(0.025 * num) > 0) {
			return ((Sample)lowSamples.get( (int)(0.025 * num)-1 )).prop;
		} 
		return 0.0;
		
	}
	
	double high95CI() {
		
		purge();
		if (keep >= (int)(0.025 * num) && (int)(0.025 * num)>0) {
			return ((Sample)hiSamples.get( (int)(0.025 * num)-1 )).prop;
		}
		return 1.0;
	
	}
	
	double relativeRepresentation() {
		
		double change = observed / (average()+epsilon);
		return (change - 1.0)*100.0;

	}

	
	double experimentalP() {
		
		if (numAbove < numBelow) {
			// Overrepresented - negative p value
			return -(numAbove+1.0)/num;
		} else {
			// Underrepresented - positive p value
			return (numBelow+1.0)/num;
		}		
	}
	
	public String getAnnotation() {
		return annotation;
	}
	
	public double getObserved() {
		return observed;
	}

    public void purge() {

	    if (hiSamples.size() > keep || lowSamples.size() > keepLow) {
    	
	        Comparator compHi = new CompareHiSample();
	        Comparator compLo = new CompareLowSample();
		    Collections.sort( hiSamples, compHi );
		    Collections.sort( lowSamples, compLo );
		    for (int i = hiSamples.size()-1; i>=keep; --i) {
		    	hiSamples.remove(i);	
		    }
		    if (hiSamples.size()>0) {
		    	lowesthigh = ((Sample)hiSamples.get(hiSamples.size()-1)).prop;
		    }
		    for (int i = lowSamples.size()-1; i>=keepLow; --i) {
		    	lowSamples.remove(i);
		    }
		    if (lowSamples.size()>0) {
		    	highestlow = ((Sample)lowSamples.get(lowSamples.size()-1)).prop;
		    }
	    }
	}

    public int isSampleSignificant( int id, int level ) {

	    // Returns +1 if sample 'id' is deemed significantly overrepresented by chance, at p value 'level'/N, 
    	// where N is the total number of samples, -1 if underrepresented, and 0 otherwise.

	    purge();
	    

	    for (int i=0; i < level; i++) {

	    	Sample s;
	    	if (hiSamples.size() > i) {
	    		s = (Sample)hiSamples.get(i);
	    		if (s.id == id) {
	    			return 1;
	    		}
	    	}
	    	if (lowSamples.size() > i) {
	    		s = (Sample)lowSamples.get(i);
	    		if (s.id == id) {
	    			return -1;
	    		}
	    	}
	    }
	    return 0;
	}


    /*
    public void dumpSamples() {

	System.out.print(getAnnotation() + "\t" );

	DecimalFormat df = new DecimalFormat();
	df.applyPattern("#.0000000");

	for (int i = 0; i < lowSamples.size(); i++) {
	    System.out.print( df.format(((Sample)(lowSamples.get(i))).prop) + "\t" );
	}

	System.out.println();

    }
    */
	
	
	private String annotation;  // Annotation identifier
	private double observed;	// observed proportion of this annotation
	private int numAbove;		// number of samples above observed proportion
	private int numBelow;		// number of samples below observed proportion
	private int num;			// total number of samples
	private double sum;			// aggregate of sampled proportions
	private double sumsq;		// aggregate of square of sampled proportions

    private List hiSamples, lowSamples;    // the 'keep' highest and lowest samples
    private double highestlow, lowesthigh;
	
}

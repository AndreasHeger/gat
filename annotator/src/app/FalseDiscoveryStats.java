package app;


import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.text.DecimalFormat;


public class FalseDiscoveryStats {
	
	
	List FDN;      // Numbers of annotations falsely designated significant over or underrep
	List FDNOver;  // Numbers of annotations falsely designated significanly overrepresented
	List FDNUnder; // Same, underrepresented
	
	int numSignificant;
	int numOverrepresented;
	int numUnderrepresented;


	private void calculateNumberSignificant( Map annotationStats, double pValueCutoff, double zThreshold ) {

		numSignificant = 0;
		numOverrepresented = 0;
		numUnderrepresented = 0;
		
		Iterator iter = annotationStats.values().iterator();
		while (iter.hasNext()) {
			Statistics stats = (Statistics)iter.next();
		
			double expPval = stats.experimentalP();
			double z = stats.z();
	    
			// criterion
			if ((Math.abs(z)>zThreshold)) {
				if (Math.abs(expPval)<=pValueCutoff) {
					if (expPval < 0.0) {
						numOverrepresented += 1;
					} else {
						numUnderrepresented += 1;
					}
					numSignificant += 1;
				}
			}
		}
    }


	
	
	private void calculateFDN( Map annotationStats, int level, int numSamples ) {
		// Calculates the false discovery number distribution, for a p-value cutoff of level/numSamples

    	FDN = new ArrayList();
    	FDNOver = new ArrayList();
    	FDNUnder = new ArrayList();
	
    	for (int id=0; id < numSamples; id++) {

    		int totFD = 0;
    		int totFDover = 0;
    		int totFDunder = 0;
    		Iterator iter = annotationStats.values().iterator();
    		while (iter.hasNext()) {
    			Statistics stats = (Statistics)iter.next();
    			int signif = stats.isSampleSignificant( id, level ); 
    			totFD += Math.abs(signif);
    			if (signif == 1) {
    				totFDover += 1;
    			}
    			if (signif == -1) {
    				totFDunder += 1;
    			}
    		}
	    
    		FDN.add( new Integer(totFD) );
    		FDNOver.add( new Integer(totFDover ));
    		FDNUnder.add( new Integer(totFDunder ));

    	}
    	
    	Collections.sort( FDN );
    	Collections.sort( FDNOver );
    	Collections.sort( FDNUnder );

	}

	


    private double avg( List l ) {

    	double tot = 0;
    	for (int i=0; i<l.size(); i++) {
    		tot += ((Integer)l.get(i)).intValue();
    	}
    	
    	return tot / l.size();
    }


    private double med( List l ) {
    	
    	return ((Integer)l.get( l.size()/2 )).intValue();

    }


    private double fiveUp( List l ) {

    	return ((Integer)l.get( (19*l.size())/20 )).intValue();
    	
    }

	
    private void printSummaryLine( int level, String category, int significant, List FDN ) {
    	
		DecimalFormat df = new DecimalFormat();

		int keep = Statistics.getKeep();
    	if (level == 0) {
    		System.out.println("--  " + category + "  (insufficiently many samples calculated)");
    	} else if (level > keep) {
	    System.out.println("--  " + category + "  " + significant + "  (required data not kept)");
	} else {
        	df.applyPattern("####0");
    		String signif = df.format( significant );
    		String med = df.format( med( FDN ));
    		String conflim = df.format( fiveUp( FDN ));
    		df.applyPattern("##0.00");
    		String avg = df.format( avg( FDN ) );
    		System.out.println ("--  " + category + "  \t" + signif + "\t\t" + avg + "\t" + med + "\t" + conflim);
    		
    	}
    	
    }


    void dumpFDRSummary( Map annotationStats, double pValueCutoff, int numSamples, double zThreshold ) {

		int level = (int)(pValueCutoff * numSamples);

    	calculateFDN( annotationStats, level, numSamples );
    	calculateNumberSignificant( annotationStats, pValueCutoff, zThreshold );

    	System.out.println("-- False Discovery summary for p-value "+pValueCutoff+":");
    	System.out.println  ("--  Category     \tObserved\tAverage\tMedian\t95% confidence limit");
    	printSummaryLine( level, "Significant", numSignificant, FDN);
    	printSummaryLine( level, "Overrep.   ", numOverrepresented, FDNOver);
    	printSummaryLine( level, "Underrep.  ", numUnderrepresented, FDNUnder);
    	
   }    
}

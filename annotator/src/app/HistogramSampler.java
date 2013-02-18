package app;


import java.util.Arrays;
import java.util.Random;


class HistogramSampler {
	
	
	public HistogramSampler( int[] histogram, int bucketSize ) {
		
		this.cdf = new double[histogram.length];
		this.bucketSize = bucketSize;
		double previous = 0;
		for (int i=0; i<histogram.length; i++) {
			cdf[i] = histogram[i] + previous;
			previous = cdf[i];
		}
		for (int i=0; i<histogram.length; i++) {
			cdf[i] /= cdf[histogram.length-1];
		}
		
	}

	
	public int sample(Random random) {
		
		double r = random.nextDouble();
		int ip = Arrays.binarySearch( cdf, r );
		int base;
		if (ip>0) {
			base = (ip-1)*bucketSize + 1;
		} else {
			base = (-(ip+1)-1)*bucketSize + 1;
		}
		if (base < 0) {
			return 0;
		}
		if (bucketSize>1) {
			return base + random.nextInt(bucketSize);
		}
		return base;		
	}
	
	
	private double[] cdf;
	private int bucketSize;
	
}
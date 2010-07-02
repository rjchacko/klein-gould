package rachele.util;

public class MathTools {
	static public double dot(double [] v1, double [] v2){
		double dotProduct = 0;
		for(int i = 0; i < v1.length; i ++)
			dotProduct += v1[i]*v2[i];
		return dotProduct;
	}

	static public double [] normalize(double [] v){
		double sum = 0;
		int L = v.length;
		double [] ret = new double [L];
		for(int i = 0; i < L; i ++)
			sum += v[i]*v[i];
		for(int i = 0; i < L; i ++)
			ret[i] = v[i]/Math.sqrt(sum);
		return ret;
	}
	
	static public double [][] normalizeRows(double [][] M){
		int L = M[0].length;
		double [][] ret = new double [L][L];
		for (int i = 0; i < L; i++)
			ret[i] = normalize(M[i]);
		return ret;
	}
	
	/**
	 * This calculates the sample variance:
	 * s^2=(1/n)*(sum(x_i-x_ave)**2)
	 * (Note that this is without Bessel's correction, which would replace n with n-1)
	 */
	static public double variance(double [] a){
		double ave = mean(a);
		double sum = 0;
		for (int i = 0; i < a.length; i++) 
			sum += Math.pow(a[i]-ave,2);
		double var = sum/(double)a.length;
		return var;
	}
	
	static public double mean(double [] a){
		double sum=0.0;
		int l = a.length;
		for (int i = 0; i < l; i++) sum += a[i];
		double average = sum/(double)l;
		return average;
	}
}

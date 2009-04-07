package chris.util;

public class MathUtil {
	
	public static double sqr(double x){
		
		return x*x;
	}
	
	public static double pow(double x, int n){
		
		double ret = x;
		
		for (int jj = 0 ; jj < n ; jj++){
			ret = x*ret;
		}
		
		return ret;
	}
	

}

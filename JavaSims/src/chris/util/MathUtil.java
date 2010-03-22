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
	
	public static int bool2bin(boolean bool){
		
		return (bool) ? 1 : 0;
	}
	
	public static int[] bool2bin(boolean bool[]){
		if (bool == null) return null;
		
		int l = bool.length;
		int[] ret = new int[l];
		
		for (int jj = 0 ; jj < l ; jj++){
			ret[jj] = (bool[jj]) ? 1 : 0;
		}
		
		return ret;
	}
	
	public static boolean bin2bool(int bin){
		
		return (bin == 1);
	}
	
	public static double[] linspace(double xmin, double xmax, int Npts){
		
		double[] ret = new double[Npts];
		for (int jj = 0 ; jj < Npts ; jj++){
			ret[jj] = xmin + jj*(xmax-xmin)/(Npts-1.);
		}
		return ret;
	}
	
	public static double[] logspace(double pmin, double pmax, int Npts){
		
		double[] ret = new double[Npts];
		for (int jj = 0 ; jj < Npts ; jj++){
			ret[jj] = Math.pow(10,pmin + jj*(pmax-pmin)/(Npts-1.));
		}
		return ret;
	}
	
	public static double[] lnspace(double pmin, double pmax, int Npts){
		
		double[] ret = new double[Npts];
		for (int jj = 0 ; jj < Npts ; jj++){
			ret[jj] = Math.exp(pmin + jj*(pmax-pmin)/(Npts-1.));
		}
		return ret;
	}
}

package chris.util;

public class CopyUtil {

	public static int[] copyArray(int[] array, int length){
		
		if(length == 0) return null;
		
		int[] ret = new int[length];
		
		for (int jj = 0 ; jj < length; jj++){
			ret[jj] = array[jj];
		}
		
		return ret;
		
	}
	
	public static double[] copyArray(double[] array, int length){
		
		if(length == 0) return null;
		
		double[] ret = new double[length];
		
		for (int jj = 0 ; jj < length; jj++){
			ret[jj] = array[jj];
		}
		
		return ret;
		
	}
	
	public static int[] copyArray(int array, int length){
		
		if(length == 0) return null;
		
		int[] ret = {array};
		
		return ret; 
	}
	
	public static double[] copyArray(double array, int length){
		
		if(length == 0) return null;
		
		double[] ret = {array};
		
		return ret; 
		
	}
	
	public static double[] invertAndScaleArray(double[] array, int length, double scale){
		
		if(length == 0) return null;
		
		double max = -100;
		double[] ret = new double[length];
		
		for (int jj = 0 ; jj < length; jj++){
			ret[jj] = 1/array[jj];
			if(ret[jj] > max) max = ret[jj];
		}
		for (int jj = 0 ; jj < length; jj++){
			ret[jj] = scale*ret[jj] / max;
		}
		
		return ret;
		
	}
	
}

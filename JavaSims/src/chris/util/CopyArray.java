package chris.util;

public class CopyArray {

	public static int[] copyArray(int[] array, int length){
		
		int[] ret = new int[length];
		
		for (int jj = 0 ; jj < length; jj++){
			ret[jj] = array[jj];
		}
		
		return ret;
		
	}
	
	public static double[] copyArray(double[] array, int length){
		
		double[] ret = new double[length];
		
		for (int jj = 0 ; jj < length; jj++){
			ret[jj] = array[jj];
		}
		
		return ret;
		
	}
	
	public static int[] copyArray(int array, int length){
		
		int[] ret = {array};
		
		return ret; 
	}
	
	public static double[] copyArray(double array, int length){
		
		double[] ret = {array};
		
		return ret; 
		
	}
	
	
}

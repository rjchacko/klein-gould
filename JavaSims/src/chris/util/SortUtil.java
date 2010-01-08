package chris.util;

import java.util.Arrays;

public class SortUtil {
	
	
	public static int[] S2LindexSort(double[] list){

		double sorted[];
		int kk, length;
		length = list.length;
		int[] ret = new int[length];
		
		sorted = CopyUtil.copyArray(list,length);
		Arrays.sort(sorted);
		for(int jj = 0 ; jj < length ; jj++){
			kk = 0;
			while(list[kk++]!=sorted[jj]);
			ret[jj] = kk - 1;
		}
		
		return ret;
	}
	
	public static double[] S2Lsort(double[] list){

		Arrays.sort(list);
		return list;
	}
	
	public static int[] S2Lsort(int[] list){

		Arrays.sort(list);
		return list;
	}
	
	public static int[] shuffleArray(int[] array, Random rand){
		
		double[] s = new double[array.length];
		int[] idx  = new int[array.length];
		int[] ret  = new int[array.length];
		
		for (int jj = 0 ; jj < array.length ; jj++){
			s[jj] = rand.nextDouble();
		}
		idx = S2LindexSort(s);
		for (int jj = 0 ; jj < array.length ; jj++){
			ret[jj] = array[idx[jj]];
		}
		
		return ret;
	}
	
	public static boolean[] shuffleArray(boolean[] array, Random rand){
		
		double[] s    = new double[array.length];
		int[] idx     = new int[array.length];
		boolean[] ret = new boolean[array.length];
		
		for (int jj = 0 ; jj < array.length ; jj++){
			s[jj] = rand.nextDouble();
		}
		idx = S2LindexSort(s);
		for (int jj = 0 ; jj < array.length ; jj++){
			ret[jj] = array[idx[jj]];
		}
		
		return ret;
	}

}

package chris.util;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;

import scikit.dataset.DatasetBuffer;
import scikit.dataset.Histogram;

public class PrintUtil {

	static public void deleteFile(String fileName){
		File file = new File(fileName);
		boolean success = file.delete();
		if (success)
			System.out.println("File deleted");
		else
			System.out.println("File delete failed");			
	}
	
	static public void printlnToFile(String fileName, String text){
		try{
			File file = new File(fileName);
			PrintWriter pw = new PrintWriter(new FileWriter(file, true), true);
			pw.println(text);
		} catch (IOException ex){
			ex.printStackTrace();
		}
	}
	
	static public void printlnToFile(String fileName, String text1, String text2){
		try{
			File file = new File(fileName);
			PrintWriter pw = new PrintWriter(new FileWriter(file, true), true);
			pw.println(text1 + "\t" + text2);
		} catch (IOException ex){
			ex.printStackTrace();
		}
	}
	
	static public void printlnToFile(String fileName, String text1, String text2, String text3){
		try{
			File file = new File(fileName);
			PrintWriter pw = new PrintWriter(new FileWriter(file, true), true);
			pw.println(text1 + "\t" + text2 + "\t" + text3);
		} catch (IOException ex){
			ex.printStackTrace();
		}
	}

	static public void printlnToFile(String fileName, String text1, String text2, String text3, String text4){
		try{
			File file = new File(fileName);
			PrintWriter pw = new PrintWriter(new FileWriter(file, true), true);
			pw.println(text1 + "\t" + text2 + "\t" + text3 + "\t" + text4);
		} catch (IOException ex){
			ex.printStackTrace();
		}
	}
	
	static public void printlnToFile(String fileName, String text1, String text2, String text3, String text4, String text5){
		try{
			File file = new File(fileName);
			PrintWriter pw = new PrintWriter(new FileWriter(file, true), true);
			pw.println(text1 + "\t" + text2 + "\t" + text3 + "\t" + text4 + "\t" + text5);
		} catch (IOException ex){
			ex.printStackTrace();
		}
	}
	
	static public void printlnToFile(String fileName, String text1, String text2, String text3, String text4, String text5, String text6){
		try{
			File file = new File(fileName);
			PrintWriter pw = new PrintWriter(new FileWriter(file, true), true);
			pw.println(text1 + "\t" + text2 + "\t" + text3 + "\t" + text4 + "\t" + text5 + "\t" + text6);
		} catch (IOException ex){
			ex.printStackTrace();
		}
	}
	
	static public void printlnToFile(String fileName, String text1, String text2, String text3, String text4, String text5, String text6, String text7){
		try{
			File file = new File(fileName);
			PrintWriter pw = new PrintWriter(new FileWriter(file, true), true);
			pw.println(text1 + "\t" + text2 + "\t" + text3 + "\t" + text4 + "\t" + text5 + "\t" + text6 + "\t" + text7);
		} catch (IOException ex){
			ex.printStackTrace();
		}
	}
	
	static public void printlnToFile(String fileName, String text1, String text2, String text3, String text4, String text5, String text6, String text7, String text8){
		try{
			File file = new File(fileName);
			PrintWriter pw = new PrintWriter(new FileWriter(file, true), true);
			pw.println(text1 + "\t" + text2 + "\t" + text3 + "\t" + text4 + "\t" + text5 + "\t" + text6 + "\t" + text7 + "\t" + text8);
		} catch (IOException ex){
			ex.printStackTrace();
		}
	}
	
	static public void printlnToFile(String fileName, String text, double data){
		try{
			File file = new File(fileName);
			PrintWriter pw = new PrintWriter(new FileWriter(file, true), true);
			pw.println(text + "\t" + data);
		} catch (IOException ex){
			ex.printStackTrace();
		}
	}
	
	static public void printlnToFile(String fileName, double d1){
		try{
			File file = new File(fileName);
			PrintWriter pw = new PrintWriter(new FileWriter(file, true), true);
			pw.println(d1);
		} catch (IOException ex){
			ex.printStackTrace();
		}
	}

	static public void printlnToFile(String fileName, double d1, double d2){
		try{
			File file = new File(fileName);
			PrintWriter pw = new PrintWriter(new FileWriter(file, true), true);
			pw.println(d1 + "\t" + d2);
		} catch (IOException ex){
			ex.printStackTrace();
		}
	}
	
	static public void printlnToFile(String fileName, String s1, double d1, double d2){
		try{
			File file = new File(fileName);
			PrintWriter pw = new PrintWriter(new FileWriter(file, true), true);
			pw.println(s1 + "\t" + d1 + "\t" + d2);
		} catch (IOException ex){
			ex.printStackTrace();
		}
	}

	static public void printlnToFile(String fileName, double d1, double d2, double d3){
		try{
			File file = new File(fileName);
			PrintWriter pw = new PrintWriter(new FileWriter(file, true), true);
			pw.println(d1 + "\t" + d2 + "\t" + d3);
		} catch (IOException ex){
			ex.printStackTrace();
		}
	}
	
	static public void printlnToFile(String fileName, double d1, double d2, double d3, double d4){
		try{
			File file = new File(fileName);
			PrintWriter pw = new PrintWriter(new FileWriter(file, true), true);
			pw.println(d1 + "\t" + d2 + "\t" + d3 + "\t" + d4);
		} catch (IOException ex){
			ex.printStackTrace();
		}
	}
	
	static public void printlnToFile(String fileName, double d1, double d2, double d3, double d4, double d5){
		try{
			File file = new File(fileName);
			PrintWriter pw = new PrintWriter(new FileWriter(file, true), true);
			pw.println(d1 + "\t" + d2 + "\t" + d3 + "\t" + d4 + "\t" + d5);
		} catch (IOException ex){
			ex.printStackTrace();
		}
	}
	
	static public void printlnToFile(String fileName, String s1, double d2, double d3, double d4, double d5){
		try{
			File file = new File(fileName);
			PrintWriter pw = new PrintWriter(new FileWriter(file, true), true);
			pw.println(s1 + "\t" + d2 + "\t" + d3 + "\t" + d4 + "\t" + d5);
		} catch (IOException ex){
			ex.printStackTrace();
		}
	}

	static public void printlnToFile(String fileName, String s1, double d2, double d3, double d4){
		try{
			File file = new File(fileName);
			PrintWriter pw = new PrintWriter(new FileWriter(file, true), true);
			pw.println(s1 + "\t" + d2 + "\t" + d3 + "\t" + d4);
		} catch (IOException ex){
			ex.printStackTrace();
		}
	}
	
	static public void printlnToFile(String fileName, double d1, double d2, double d3, double d4, double d5, double d6){
		try{
			File file = new File(fileName);
			PrintWriter pw = new PrintWriter(new FileWriter(file, true), true);
			pw.println(d1 + "\t" + d2 + "\t" + d3 + "\t" + d4 + "\t" + d5 + "\t" + d6);
		} catch (IOException ex){
			ex.printStackTrace();
		}
	}
	
	static public void printlnToFile(String fileName, double d1, double d2, double d3, double d4, double d5, double d6, double d7){
		try{
			File file = new File(fileName);
			PrintWriter pw = new PrintWriter(new FileWriter(file, true), true);
			pw.println(d1 + "\t" + d2 + "\t" + d3 + "\t" + d4 + "\t" + d5 + "\t" + d6 + "\t" + d7);
		} catch (IOException ex){
			ex.printStackTrace();
		}
	}

	static public void printlnToFile(String fileName, double d1, double d2, double d3, double d4, double d5, double d6, double d7, double d8){
		try{
			File file = new File(fileName);
			PrintWriter pw = new PrintWriter(new FileWriter(file, true), true);
			pw.println(d1 + "\t" + d2 + "\t" + d3 + "\t" + d4 + "\t" + d5 + "\t" + d6 + "\t" + d7 + "\t" + d8);
		} catch (IOException ex){
			ex.printStackTrace();
		}
	}
	
	static public void printWalkerData(String fileName, int[] array, int N, double offSet, double binWidth){
		try{
			File file = new File(fileName);
			PrintWriter pw = new PrintWriter(new FileWriter(file, true), true);
			for (int ii = 0 ; ii < N ; ii++){				
				pw.print(ii);
				pw.print("\t");
				pw.print(offSet + ii*binWidth);
				pw.print("\t");
				pw.print(array[ii]);
				pw.print("\t");
				pw.print((double)(array[ii])/(offSet + ii*binWidth));
				pw.println();
			}
			pw.close();
		}
		catch (IOException ex){
			ex.printStackTrace();
		}
	}
	
	static public void printWalkerData(String fileName, double[] array, int N, double offSet, double binWidth){
		try{
			File file = new File(fileName);
			PrintWriter pw = new PrintWriter(new FileWriter(file, true), true);
			for (int ii = 0 ; ii < N ; ii++){				
				pw.print(ii);
				pw.print("\t");
				pw.print(offSet + ii*binWidth);
				pw.print("\t");
				pw.print(array[ii]);
				pw.print("\t");
				pw.print((double)(array[ii])/(offSet + ii*binWidth));
				pw.println();
			}
			pw.close();
		}
		catch (IOException ex){
			ex.printStackTrace();
		}
	}

	
	static public void printWalkerData(String fileName, int[][] array, int N){
		try{
			File file = new File(fileName);
			PrintWriter pw = new PrintWriter(new FileWriter(file, true), true);
			for (int ii = 0 ; ii < N ; ii++){				
				for(int jj = 0 ; jj < N ; jj++){
					pw.print((double)(ii)/(double)(N)-0.5);
					pw.print("\t");
					pw.print((double)(jj)/(double)(N)-0.5);
					pw.print("\t");
					pw.print(array[ii][jj]);
					pw.println();
				}
				pw.println();
			}
			pw.close();
		}
		catch (IOException ex){
			ex.printStackTrace();
		}		
	}
	
	
	static public void printVectorAsMatrixToFile(String fileName, double[] array, int m, int n){
		try{
			
			File file = new File(fileName);
			PrintWriter pw = new PrintWriter(new FileWriter(file, true), true);

			for (int ii = 0 ; ii < m*n ; ii++){
				pw.print(array[ii]);
				if ( (ii + 1)%n == 0){
					pw.println();
				}
				else{
					pw.print("\t");
				}
			}
			
			pw.close();

		} catch (IOException ex){
			ex.printStackTrace();
		}
	}
	
	static public void printVectorAsMatrixToFile(String fileName, int[] array, int m, int n){
		try{
			
			File file = new File(fileName);
			PrintWriter pw = new PrintWriter(new FileWriter(file, true), true);

			for (int ii = 0 ; ii < m*n ; ii++){
				pw.print(array[ii]);
				if ( (ii + 1)%n == 0){
					pw.println();
				}
				else{
					pw.print("\t");
				}
			}
			
			pw.close();

		} catch (IOException ex){
			ex.printStackTrace();
		}
	}
	
	
	static public void printMxNarrayToFile(String fileName, double[][] array, int m, int n){
		try{
			
			File file = new File(fileName);
			PrintWriter pw = new PrintWriter(new FileWriter(file, true), true);
	
			for (int ii = 0 ; ii < m ; ii++){
				for(int jj = 0 ; jj < n ; jj++){
					pw.print(array[jj][ii]);
					pw.print("\t");
				}
				pw.println();
			}
			
			pw.close();

		} catch (IOException ex){
			ex.printStackTrace();
		}
	}
	
	static public void printArray4Splot(String fileName, double[][] array, int Yindx, int Zindx, int N, double Xval){
		try{
			File file = new File(fileName);
			PrintWriter pw = new PrintWriter(new FileWriter(file, true), true);
			
		
		for (int ii = 0 ; ii < N ; ii++){
			pw.print(Xval);
			pw.print("\t");
			pw.print(array[Yindx][ii]);
			pw.print("\t");
			pw.print(array[Zindx][ii]);
			pw.println();
		}
		pw.println();
		
		pw.close();
		}
		catch (IOException ex){
		ex.printStackTrace();
	}
		
	}
	
	static public void printScalarsAndVectorToFile(String fileName, double t1, double t2, double[] vector){

		int Length = vector.length;	
		try{
			File file = new File(fileName);
			PrintWriter pw = new PrintWriter(new FileWriter(file, true), true);
			pw.print(t1);
			pw.print("\t");
			pw.print(t2);
			pw.print("\t");
			for (int ii = 0 ; ii<Length ; ii++){
				pw.print(vector[ii]);
				pw.print("\t");
			}
			pw.println();

		} catch (IOException ex){
			ex.printStackTrace();
		}
	}
	
	static public void printScalarAndVectorToFile(String fileName, double t1, double[] vector){

		int Length = vector.length;	
		try{
			File file = new File(fileName);
			PrintWriter pw = new PrintWriter(new FileWriter(file, true), true);
			pw.print(t1);
			pw.print("\t");
			for (int ii = 0 ; ii<Length ; ii++){
				pw.print(vector[ii]);
				pw.print("\t");
			}
			pw.println();

		} catch (IOException ex){
			ex.printStackTrace();
		}
	}
	
	static public void printScalarAndVectorToFile(String fileName, double t1, int[] vector){

		int Length = vector.length;	
		try{
			File file = new File(fileName);
			PrintWriter pw = new PrintWriter(new FileWriter(file, true), true);
			pw.print(t1);
			pw.print("\t");
			for (int ii = 0 ; ii<Length ; ii++){
				pw.print(vector[ii]);
				pw.print("\t");
			}
			pw.println();

		} catch (IOException ex){
			ex.printStackTrace();
		}
	}
	

	static public void printlnToFile(String fileName, Boolean b1, Boolean b2){
		try{
			File file = new File(fileName);
			PrintWriter pw = new PrintWriter(new FileWriter(file, true), true);
			if(b1 = true){
				pw.print("true");
				pw.print("\t");
			}
			else{
				pw.print("false");
				pw.print("\t");
			}
			if(b2 = true){
				pw.println("true");
			}
			else{
				pw.println("false");
			}
			
			
		} catch (IOException ex){
			ex.printStackTrace();
		}
	}

	static public void printHistToFile(String fileName, Histogram h){
		DatasetBuffer hh = h.copyData();
		try{
			File file = new File(fileName);
			PrintWriter pw = new PrintWriter(new FileWriter(file, true), true);
			for (int jj = 0 ; jj < hh.size() ; jj++){
				pw.print(hh.x(jj));
				pw.print("\t");
				pw.println(hh.y(jj));
			}
		} catch (IOException ex){
			ex.printStackTrace();
		}
	}

	public static void printArraysToFile(String fileName, int[] a1, int[] a2) {
		int Length = a1.length;	
		if(a2.length != Length) throw new IllegalArgumentException("Arrays have different lengths");
		try{
			File file = new File(fileName);
			PrintWriter pw = new PrintWriter(new FileWriter(file, true), true);
			for (int ii = 0 ; ii<Length ; ii++){
				pw.print(a1[ii]);
				pw.print("\t");
				pw.println(a2[ii]);
			}
		} catch (IOException ex){
			ex.printStackTrace();
		}
	}
	
	public static void printArraysToFile(String fileName, double[] a1, double[] a2) {
		int Length = a1.length;	
		if(a2.length != Length) throw new IllegalArgumentException("Arrays have different lengths");
		try{
			File file = new File(fileName);
			PrintWriter pw = new PrintWriter(new FileWriter(file, true), true);
			for (int ii = 0 ; ii<Length ; ii++){
				pw.print(a1[ii]);
				pw.print("\t");
				pw.println(a2[ii]);
			}
		} catch (IOException ex){
			ex.printStackTrace();
		}
	}

	public static void printVectorToFile(String fileName, double[] vec, int Length){
		try{
			File file = new File(fileName);
			PrintWriter pw = new PrintWriter(new FileWriter(file, true), true);
			for (int ii = 0 ; ii<Length ; ii++){
				pw.println(vec[ii]);
			}
		} catch (IOException ex){
			ex.printStackTrace();
		}
		
	}
	
	public static void printVectorToFile(String fileName, int[] vec, int Length){
		try{
			File file = new File(fileName);
			PrintWriter pw = new PrintWriter(new FileWriter(file, true), true);
			for (int ii = 0 ; ii<Length ; ii++){
				pw.println(vec[ii]);
			}
		} catch (IOException ex){
			ex.printStackTrace();
		}
		
	}
	
	public static void printVectorToFile(String fileName, int[] vec){
		int Length = vec.length;
		try{
			File file = new File(fileName);
			PrintWriter pw = new PrintWriter(new FileWriter(file, true), true);
			for (int ii = 0 ; ii<Length ; ii++){
				pw.println(vec[ii]);
			}
		} catch (IOException ex){
			ex.printStackTrace();
		}
		
	}
	
	public static void printVectorToFile(String fileName, double[] vec){
		int Length = vec.length;
		try{
			File file = new File(fileName);
			PrintWriter pw = new PrintWriter(new FileWriter(file, true), true);
			for (int ii = 0 ; ii<Length ; ii++){
				pw.println(vec[ii]);
			}
		} catch (IOException ex){
			ex.printStackTrace();
		}
		
	}
	
	
}

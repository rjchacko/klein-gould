package chris.util;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;

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
	
	static public void printArrayToFile(String fileName, int[] array, int m, int n){
		try{
			
			File file = new File(fileName);
			PrintWriter pw = new PrintWriter(new FileWriter(file, true), true);
			
			
			// print info here
//			pw.println("m by n array");
//			pw.println("printing:");
//			pw.println("i = 0, 1, 2, . . . m-1 ");
//			pw.println("i = m, m+1 , . . . ");
//			pw.println("    .");
//			pw.println("    .");
//			pw.println("    .");
//			pw.println("    n, n+1, . . . n+m");
//			pw.print("m = ");
//			pw.println(m);
//			pw.print("n = ");
//			pw.println(n);
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
	
	static public void printArrayToFile(String fileName, double[] array, int m, int n){
		try{
			
			File file = new File(fileName);
			PrintWriter pw = new PrintWriter(new FileWriter(file, true), true);
			
			
			// print info here
			pw.println("m by n array");
			pw.println("printing:");
			pw.println("i = 0, 1, 2, . . . m-1 ");
			pw.println("i = m, m+1 , . . . ");
			pw.println("    .");
			pw.println("    .");
			pw.println("    .");
			pw.println("    n, n+1, . . . n+m");
			pw.print("m = ");
			pw.println(m);
			pw.print("n = ");
			pw.println(n);
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
	
	static public void printTimeAndVectorToFile(String fileName, double T, double[] vector){

		int Length = vector.length;	
		try{
			File file = new File(fileName);
			PrintWriter pw = new PrintWriter(new FileWriter(file, true), true);
			pw.print(T);
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
	
	/*
	 * 
	 * Replace this with the above method!!!!!
	 * 
	 */
	
	static public void print2TimeAndVectorToFile(String fileName, double T1, double T2, double[] vector){

		int Length = vector.length;	
		try{
			File file = new File(fileName);
			PrintWriter pw = new PrintWriter(new FileWriter(file, true), true);
			pw.print(T1);
			pw.print("\t");
			pw.print(T2);
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
	
	static public void print2TimeAndVectorToFile(String fileName, double T1, double T2, int[] vector){

		int Length = vector.length;	
		try{
			File file = new File(fileName);
			PrintWriter pw = new PrintWriter(new FileWriter(file, true), true);
			pw.print(T1);
			pw.print("\t");
			pw.print(T2);
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
	
	// this is ad hoc...delete it!!!!
	static public void printTimeAndVectorToFile(String fileName, double t1, int[][] array, int index, int length){
		int Length = length;	
		try{
			File file = new File(fileName);
			PrintWriter pw = new PrintWriter(new FileWriter(file, true), true);
			pw.print(t1);
			pw.print("\t");
			for (int ii = 0 ; ii<Length ; ii++){
				pw.print(array[index][ii]);
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
	
	
}

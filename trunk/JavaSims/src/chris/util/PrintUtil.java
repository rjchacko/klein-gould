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
	
	static public void printlnToFile(String fileName, String text, double data){
		try{
			File file = new File(fileName);
			PrintWriter pw = new PrintWriter(new FileWriter(file, true), true);
			pw.println(text + "\t" + data);
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
			pw.println(d1 + "\t" + d2 + "\t" + d3 + "\t" + d4 + "\t" + d5+ "\t" + d6);
		} catch (IOException ex){
			ex.printStackTrace();
		}
	}

}

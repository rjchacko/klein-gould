package rachele.util;

import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.EOFException;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;

import scikit.dataset.Accumulator;
import scikit.dataset.DatasetBuffer;
import scikit.jobs.params.Parameters;

public class FileUtil {
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

	static public void printlnToFile(String fileName, String text, double data){
		try{
			File file = new File(fileName);
			PrintWriter pw = new PrintWriter(new FileWriter(file, true), true);
			pw.println(text + " " + data);
		} catch (IOException ex){
			ex.printStackTrace();
		}
	}

	static public void printlnToFile(String fileName, String text, String text2){
		try{
			File file = new File(fileName);
			PrintWriter pw = new PrintWriter(new FileWriter(file, true), true);
			pw.println(text + " " + text2);
		} catch (IOException ex){
			ex.printStackTrace();
		}
	}
	
	static public void printlnToFile(String fileName, double d1, double d2){
		try{
			File file = new File(fileName);
			PrintWriter pw = new PrintWriter(new FileWriter(file, true), true);
			pw.println(d1 + " " + d2);
		} catch (IOException ex){
			ex.printStackTrace();
		}
	}

	static public void printlnToFile(String fileName, double d1, double d2, double d3){
		try{
			File file = new File(fileName);
			PrintWriter pw = new PrintWriter(new FileWriter(file, true), true);
			pw.println(d1 + " " + d2 + " " + d3);
		} catch (IOException ex){
			ex.printStackTrace();
		}
	}

	static public void printlnToFile(String fileName, double d1, double d2, double d3, double d4){
		try{
			File file = new File(fileName);
			PrintWriter pw = new PrintWriter(new FileWriter(file, true), true);
			pw.println(d1 + " " + d2 + " " + d3 + " " + d4);
		} catch (IOException ex){
			ex.printStackTrace();
		}
	}

	static public void printAccumToFile(String fileName, Accumulator acc){
		DatasetBuffer data = acc.copyData();
		int size = data.size();
		try{
			File file = new File(fileName);
			PrintWriter pw = new PrintWriter(new FileWriter(file, true), true);
			for (int i = 0; i < size; i++)
				pw.println(data.x(i) + " " + data.y(i) + " " + data.errorY(i));
		} catch (IOException ex){
			ex.printStackTrace();
		}
	}
	
	static public void writeConfigToFile(String FileName, int size, double [] A){
	// write configuration A[] in binary
		try {
			File pathFile = new File(FileName);
			DataOutputStream dos = new DataOutputStream(new FileOutputStream(pathFile, true));
			for (int i = 0; i < size; i ++){
				dos.writeInt(i);
				dos.writeChar('\t');
				dos.writeDouble(A[i]);
				//System.out.println(A[i]);
				dos.writeChar('\n');
			}
			dos.close();
		} catch (IOException ex){
			ex.printStackTrace();
		}
		System.out.println("config written");
	}

	static public double [] readConfigFromFile(String FileName, int size){
		
		double [] A = new double [size];
		try{
			File myFile = new File(FileName);
			DataInputStream dis = new DataInputStream(new FileInputStream(myFile));
			try{
				while(true){
					int i = dis.readInt();
					dis.readChar();       // throws out the tab
					A[i] = dis.readDouble();
					dis.readChar();
				}
			} catch (EOFException e) {
			}

		} catch (FileNotFoundException e) {
			System.err.println("FileStreamsTest: " + e);
		} catch (Exception ex) {
			ex.printStackTrace();
		}	
		
		return A;
	}
	
	static public void initFile(String fileName, Parameters params){
		deleteFile(fileName);
		try{
			File file = new File(fileName);
			PrintWriter pw = new PrintWriter(new FileWriter(file, true), true);
			String [] keys = params.keys();
			for (int i = 0; i < keys.length; i++)
				pw.println("# " + keys[i] + " " + params.sget(keys[i]));
		} catch (IOException ex){
			ex.printStackTrace();
		}
	}

	static public void initFile(String fileName, Parameters params, String message1){
		deleteFile(fileName);
		try{
			File file = new File(fileName);
			PrintWriter pw = new PrintWriter(new FileWriter(file, true), true);
			pw.println(message1);
			String [] keys = params.keys();
			for (int i = 0; i < keys.length; i++)
				pw.println("# " + keys[i] + " " + params.sget(keys[i]));
		} catch (IOException ex){
			ex.printStackTrace();
		}
	}
	
	static public void initFile(String fileName, Parameters params, String message1, String message2){
		deleteFile(fileName);
		try{
			File file = new File(fileName);
			PrintWriter pw = new PrintWriter(new FileWriter(file, true), true);
			pw.println(message1);
			pw.println(message2);
			String [] keys = params.keys();
			for (int i = 0; i < keys.length; i++)
				pw.println(" # " + keys[i] + " " + params.sget(keys[i]));
		} catch (IOException ex){
			ex.printStackTrace();
		}
	}
	
}

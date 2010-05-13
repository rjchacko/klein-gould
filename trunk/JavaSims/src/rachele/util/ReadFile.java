package rachele.util;

import java.io.BufferedInputStream;
import java.io.BufferedReader;
import java.io.DataInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;


public class ReadFile {
	/**
	 * This method reads a text file line by line and print to the console. It uses
	 * FileOutputStream to read the file.
	 * 
	 */
	@SuppressWarnings("deprecation")
	public static void printFile(String filename){
		File file = new File(filename);
		FileInputStream fis = null;
		BufferedInputStream bis = null;
		DataInputStream dis = null;

		try {
			fis = new FileInputStream(file);

			// Here BufferedInputStream is added for fast reading.
			bis = new BufferedInputStream(fis);
			dis = new DataInputStream(bis);

			// dis.available() returns 0 if the file does not have more lines.
			while (dis.available() != 0) {

				// this statement reads the line from the file and print it to
				// the console.
				String readLine = dis.readLine();
				System.out.println(readLine);
			}

			// dispose all the resources after using them.
			fis.close();
			bis.close();
			dis.close();

		} catch (FileNotFoundException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	/**
	 * This method reads a text file line by line and returns
	 * the specified line as a string.
	 */
	@SuppressWarnings("deprecation")
	public static String readWord(String filename, int getLine, int getWord){
		File file = new File(filename);
		FileInputStream fis = null;
		BufferedInputStream bis = null;
		DataInputStream dis = null;
		String line = "0";
		String[] words;

		try {
			fis = new FileInputStream(file);

			// Here BufferedInputStream is added for fast reading.
			bis = new BufferedInputStream(fis);
			dis = new DataInputStream(bis);

			
			int lineNo = 0;
			// dis.available() returns 0 if the file does not have more lines.
			while (dis.available() != 0) {
				lineNo += 1;
				// this statement reads the line from the file and print it to
				// the console.
				String readLine = dis.readLine();
				if (lineNo == getLine) line = readLine;
			}

			// dispose all the resources after using them.
			fis.close();
			bis.close();
			dis.close();

		} catch (FileNotFoundException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		}
		words = line.split("\\s+");
		String word = words[getWord];
		return word;
	}
	
	/**
	 * This method reads a text file line by line and returns
	 * the specified line as a string.
	 */
	@SuppressWarnings("deprecation")
	public static String readLine(String filename, int lineToGet){
		File file = new File(filename);
		FileInputStream fis = null;
		BufferedInputStream bis = null;
		DataInputStream dis = null;
		String ret = "0";

		try {
			fis = new FileInputStream(file);

			// Here BufferedInputStream is added for fast reading.
			bis = new BufferedInputStream(fis);
			dis = new DataInputStream(bis);

			
			int line = 0;
			// dis.available() returns 0 if the file does not have more lines.
			while (dis.available() != 0) {
				line += 1;
				// this statement reads the line from the file and print it to
				// the console.
				String readLine = dis.readLine();
				if (line == lineToGet) ret = readLine;
			}

			// dispose all the resources after using them.
			fis.close();
			bis.close();
			dis.close();

		} catch (FileNotFoundException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		}
		return ret;
	}
	
	
	/**
	 * This method prints the whole file into one string.
	 * 
	 */
	public static String fileToString(String fileName) {
        File file = new File(fileName);
        StringBuilder contents = new StringBuilder();
        BufferedReader input = null;
        try {
            input = new BufferedReader(new FileReader(file));
            String line = null;
            while((line = input.readLine()) != null)
                contents.append(line);
            return contents.toString();
        } catch (Exception e) {
            throw new RuntimeException(e);
        } finally {
            try {
                if(input != null) input.close();
            } catch(Exception e) {
            }
        }
    }
}



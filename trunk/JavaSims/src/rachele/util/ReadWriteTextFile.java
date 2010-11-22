package rachele.util;

import java.io.*;

public class ReadWriteTextFile {

  /**
  * Fetch the entire contents of a text file, and return it in a String.
  * This style of implementation does not throw Exceptions to the caller.
  *
  * @param aFile is a file which already exists and can be read.
  */
  static public String getContents(File aFile) {
    //...checks on aFile are elided
    StringBuilder contents = new StringBuilder();
    
    try {
      //use buffering, reading one line at a time
      //FileReader always assumes default encoding is OK!
      BufferedReader input =  new BufferedReader(new FileReader(aFile));
      try {
        String line = null; //not declared within while loop
        /*
        * readLine is a bit quirky :
        * it returns the content of a line MINUS the newline.
        * it returns null only for the END of the stream.
        * it returns an empty String if two newlines appear in a row.
        */
        while (( line = input.readLine()) != null){
          contents.append(line);
          contents.append(System.getProperty("line.separator"));
        }
      }
      finally {
        input.close();
      }
    }
    catch (IOException ex){
      ex.printStackTrace();
    }
    
    return contents.toString();
  }
  
  /**
   * Fetch the entire contents of a line of a text file, and return it in a String.
   * This style of implementation does not throw Exceptions to the caller.
   *
   * @param aFile is a file which already exists and can be read.
   */
   static public String getLineContents(File aFile, int lineNo) {
     //...checks on aFile are elided
     StringBuilder contents = new StringBuilder();
     int currentline = 0;
     
     try {
       //use buffering, reading one line at a time
       //FileReader always assumes default encoding is OK!
       BufferedReader input =  new BufferedReader(new FileReader(aFile));
       try {
         String line = null; //not declared within while loop
         /*
         * readLine is a bit quirky :
         * it returns the content of a line MINUS the newline.
         * it returns null only for the END of the stream.
         * it returns an empty String if two newlines appear in a row.
         */
         while (( line = input.readLine()) != null){
        	 if(currentline == lineNo){
        		 contents.append(line);
        		 contents.append(System.getProperty("line.separator"));
        	 }
        	 currentline+=1;
         }
       }
       finally {
         input.close();
       }
     }
     catch (IOException ex){
       ex.printStackTrace();
     }
     
     return contents.toString();
   }
   
   /**
    * Fetches the last word from a line of a text file, and returns it in a String.
    * This style of implementation does not throw Exceptions to the caller.
    *
    *Appropriately removes last character at end of line
    *
    * @param aFile is a file which already exists and can be read.
    */
    static public String getLastWordOfLine(File aFile, int lineNo) {
    	
    	String lineContents = getLineContents(aFile, lineNo);
    	String[] words = lineContents.split(" ");
		int noWords = words.length;
//		String lastword = words[noWords-1];//= words[noWords-1].substring(0, words[noWords-1].length()-1);
		String lastword = words[noWords-1].substring(0, words[noWords-1].length()-1);
      return lastword;
    }

    /**
     * Fetches a requested word from a line of a text file, and returns it in a String.
     * This style of implementation does not throw Exceptions to the caller.
     *
     * (Does not remove last character at end of word.)
     *
     * @param aFile is a file which already exists and can be read.
     */
    static public String getWordOfLine(File aFile, int lineNo, int wordno) {
    	
    	String lineContents = getLineContents(aFile, lineNo);
    	String[] words = lineContents.split(" ");
      return words[wordno];
    }
    
    /**
     * Fetches the last word from a line of a text file, and returns it in a String.
     * This style of implementation does not throw Exceptions to the caller.
     *
     *Appropriately removes last character at end of line
     *
     * @param aFile is a file which already exists and can be read.
     */
    static public String getLastWordsOfLine(File aFile, int lineNo, int wordAfter) {
    	
        StringBuilder contents = new StringBuilder();
		int noWords = getLineContents(aFile, lineNo).split(" ").length;
		for(int i = wordAfter; i<noWords; i++){

			if(i!=noWords-1){
				contents.append(getWordOfLine(aFile, lineNo, i));
				contents.append(" ");
			}
			else contents.append(getLastWordOfLine(aFile, lineNo)); 
		}
		String ret = contents.toString();
      return ret;
    }
    
}
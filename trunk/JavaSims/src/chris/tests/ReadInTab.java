package chris.tests;

public class ReadInTab {
	static final public void main(final String[] args) {
		
		String test = "Hello World=Java Program=Wonder if this works";
		// Portion before a delimiter (=)
		int posDelimiter = test.indexOf("=");
		System.out.println(test.substring(0, posDelimiter));
		// Portion after a delimiter (=)
		//System.out.println(test.substring(posDelimiter + 1));
		String test2=test.substring(posDelimiter + 1);
		System.out.println(test2.substring(0, posDelimiter));
		System.out.println(test2.substring(posDelimiter+2));		
		// Portion between two delimiters (space)
		//int posStart = test.indexOf(" ");
		//int posFinish = test.indexOf(" ", posStart + 1);
		//System.out.println(test.substring(posStart + 1, posFinish));
		
		
		test = "Hello World" + "\t" + "Java Program" + "\t" +"Wonder if this works";
		// Portion before a delimiter (=)
		posDelimiter = test.indexOf("\t");
		System.out.println(test.substring(0, posDelimiter));
		// Portion after a delimiter (=)
		//System.out.println(test.substring(posDelimiter + 1));
		test2=test.substring(posDelimiter + 1);
		System.out.println(test2.substring(0, posDelimiter));
		System.out.println(test2.substring(posDelimiter+2));		
		// Portion between two delimiters (space)
		//int posStart = test.indexOf(" ");
		//int posFinish = test.indexOf(" ", posStart + 1);
		//System.out.println(test.substring(posStart + 1, posFinish));
		
	}
}
//String fin = model.outdir + File.separator+"Test2.txt";
//
//try { 
//
//    File f = new File(fin); 
//    FileInputStream fis = new FileInputStream(f); 
//    BufferedInputStream bis = new BufferedInputStream(fis); 
//    DataInputStream dis = new DataInputStream(bis);  
//
//    while ( (record=dis.readLine()) != null ) { 
//       System.out.println(record); 
//       foo[ii++]=Double.parseDouble(record);
//       // Integer.parseInt() also exists!
//    } 
//
// } catch (IOException e) { 
//    // catch io errors from FileInputStream or readLine() 
//    System.out.println("ERROR!" + e.getMessage()); 
//
// }
// 
// 
// for ( int kk = 0 ; kk<8 ; kk++){
//	 foo[kk]=Math.sqrt(foo[kk]);
// }
// 
// for ( int kk = 0 ; kk<8 ; kk++){
//	 System.out.println(foo[kk]);
// }
// 
package chris.tests;

import chris.util.DirUtil;

public class DirTest {

	public static void main(String[] args) {
		
		String testpath1 = "/Users/cserino/Desktop/TestHere/";
		String testpath2 = "/Users/cserino/Desktop/TestThere/";
		String testpath3 = "/Users/cserino/Desktop/TopDir/SubDir/";
		
		System.out.println("Testing");
		DirUtil.MkDir(testpath1);
		DirUtil.MkDir(testpath2);
		//DirUtil.MkDir(testpath3);
		DirUtil.MkDirs(testpath3);
		System.out.println("Done");
		
	}

}

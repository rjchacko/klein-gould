package chris.util;

import java.io.File;

public class DirUtil {
	
	public static void MkDir(String PATH){
		
		if (new File(PATH).exists()){
			System.out.println("Directory already exists.");
			return;
		}
		
		try {
			boolean mddir = (new File(PATH) ).mkdir();
			if (mddir) {
			      System.out.println("Directory: " + PATH + " created");
			}    
		}
		catch (Exception e) {
			System.err.println("Error: " + e.getMessage());
		}

		return;
	}
	
	public static void MkDirs(String PATH){
		
		if (new File(PATH).exists()){
			System.out.println("Directory already exists.");
			return;
		}
		
		try {
			boolean mddir = (new File(PATH)).mkdirs();
			if (mddir) {
			      System.out.println("Directory: " + PATH + " created");
			}    
		}
		catch (Exception e) {
			System.err.println("Error: " + e.getMessage());
		}

		return;
	}

}

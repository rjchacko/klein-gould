package chris.util;

import java.io.File;
import java.io.FileFilter;
import java.io.FilenameFilter;

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
	
	public static File[] getFiles(String din){
		File dir     = new File(din);
		FileFilter fileFilter = new FileFilter(){
			public boolean accept(File file) {
				return file.isFile();
			}
		};

		return dir.listFiles(fileFilter);
	}
	
	public static File[] getFiles(String din, final String fltr){
		File dir     = new File(din);
		FilenameFilter fileFilter = new FilenameFilter(){
			public boolean accept(File file, String name) {
				return name.contains(fltr);
			}
		};
		
		return dir.listFiles(fileFilter);
	}

}

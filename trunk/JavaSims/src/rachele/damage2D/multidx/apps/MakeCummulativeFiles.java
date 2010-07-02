package rachele.damage2D.multidx.apps;

import scikit.jobs.Control;
import scikit.jobs.Simulation;
import scikit.jobs.params.DirectoryValue;
import java.io.File;

public class MakeCummulativeFiles extends Simulation{

	public static void main(String[] args) {
		new Control(new MakeCummulativeFiles(), "Make Cummulative Files");
	}
	
	public void animate() {
		
	}

	public void clear() {
		
	}

	public void load(Control c) {
		params.add("Data Dir",new DirectoryValue("/Users/erdomi/data/damage/contract2/testRuns/"));
		
	}

	
	public void run() {
		String parentDir = params.sget("Data Dir");
		String shFile = parentDir + "sh.txt";
		findFiles();
	}
	
	/**
	* Find any files that begin with sh
	*/
	void findFiles(){
		    File folder = new File(params.sget("Data Dir"));
		    File[] listOfFiles = folder.listFiles();
		    int noOfFiles = listOfFiles.length;

		    for (int i = 0; i < noOfFiles; i++) {
		      if (listOfFiles[i].isFile()) {
		        System.out.println("File " + listOfFiles[i].getName());
		      } else if (listOfFiles[i].isDirectory()) {
		        System.out.println("Directory " + listOfFiles[i].getName());
		      }
		    }
		 

	}

}

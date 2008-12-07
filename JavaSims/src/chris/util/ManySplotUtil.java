package chris.util;

import scikit.dataset.Histogram;
import scikit.jobs.Control;
import scikit.jobs.Simulation;
import scikit.jobs.params.FileValue;

public class ManySplotUtil extends Simulation{

	Histogram histGLOB;
	
	public static void main(String[] args) {
		new Control(new ManySplotUtil(), "Create for GNUPLOT Splot from TAB");
	}
	
	public void load(Control c) {
		
		params.add("TAB File", new FileValue("/Users/cserino/Desktop/"));
		params.add("Number of Header Rows / Rows to Skip", (int) 0);
		params.add("Column Number for Y Data", (int) 2);
		params.add("Column Number for Z Data", (int) 3);
		params.add("Status");
				
	}
	
	public void animate() {
	}

	public void clear() {
	}

	public void run() {
		
		
		ReadInUtil rin = new ReadInUtil(params.sget("TAB File"));
		
		double[][] tmp = rin.getData(new int[] {1,2,3},params.iget("Number of Header Rows / Rows to Skip"));
		
		
		String tmpfout = "/Users/cserino/Desktop/TESTOUT.txt";
		
		
		PrintUtil.printArrayToFile(tmpfout,tmp,tmp[0].length,tmp.length);
		
		return;
	}
	

}

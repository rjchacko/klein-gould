package chris.tests;

import chris.util.PrintUtil;
import chris.util.ReadInUtil;
import scikit.jobs.Control;
import scikit.jobs.Simulation;
import scikit.jobs.params.FileValue;

public class ReadInTest extends Simulation{

	private String fin;
	private ReadInUtil imprt;
	private double ret[][];
	
	public static void main(String[] args) {
		new Control(new ReadInTest(), "Test LinFit");
	}
	
	public void animate() {

		return;
	}

	public void clear() {

		return;
	}

	public void load(Control c) {
		params.add("Input File",new FileValue("/Users/cserino/Desktop/TestLinFit.txt"));
		return;
	}

	public void run() {
			
		fin = params.sget("Input File");
		
		imprt = new ReadInUtil(fin);
		
		ret = imprt.getData(new int[]{3,1,6}, 1);
		
		for(int jj = 0 ; jj < ret[0].length ; jj++){
			PrintUtil.printlnToFile("/Users/cserino/Desktop/OutputTest.txt",ret[0][jj], ret[1][jj], ret[2][jj]);
		}
			
	}
	


}

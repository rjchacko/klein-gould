package chris.tests;

import scikit.jobs.Control;
import scikit.jobs.Simulation;
import scikit.jobs.params.FileValue;
import chris.util.FitUtil;
import chris.util.ReadInUtil;

public class ReadInTest extends Simulation{

	private String fin;
	private ReadInUtil imprt;
	@SuppressWarnings("unused")
	private double rin[][], ret[];
	@SuppressWarnings("unused")
	private FitUtil fitter;
	
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
		params.add("Input File",new FileValue("/Users/cserino/Desktop"));
		return;
	}

	public void run() {
			
		fin = params.sget("Input File");
		
		imprt = new ReadInUtil(fin);
		
		rin = imprt.getData(new int[]{1,8}, 1);
		
		double[] tempx = new double[rin[0].length];
		double[] tempy = new double[rin[0].length];

		for (int jj = 0 ; jj < rin[0].length ; jj++){
			tempx[jj] = rin[0][jj];
			tempy[jj] = rin[1][0]/rin[1][jj];
		}
		
		
//		fitter = new FitUtil(rin[1].length);
//		
//		ret = fitter.fit(tempx,tempy,30*30,false);
//		
//		System.out.print("chi2 = ");
//		System.out.println(ret[4]);
//		
//			
	}
	


}

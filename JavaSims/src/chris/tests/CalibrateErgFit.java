package chris.tests;

import scikit.jobs.Control;
import scikit.jobs.Job;
import scikit.jobs.Simulation;
import scikit.jobs.params.FileValue;
import chris.util.CopyUtil;
import chris.util.ReadInUtil;

public class CalibrateErgFit extends Simulation{

	@SuppressWarnings("unused")
	private double xvals[], yvals[], m, b, chi2;
	private ReadInUtil rin;
	//private FitUtil fitter;

	
	public static void main(String[] args) {
		new Control(new CalibrateErgFit(), "Calibrate Ergoduc Fitter");
	}
	
	public void animate() {

		return;
	}

	
	public void clear() {

		return;
	}

	
	public void load(Control c) {
		params.add("Input File (tab del.)",new FileValue("/Users/cserino/Research/Data2/ErgPhaseDiagram/NearestNeighbor/TestErg.txt"));
		params.add("Time Col. Num.", (int) 1);
		params.add("Metric Col. Num.", (int) 8);
	}

	public void run() {
		
		int count = 0;
		
		double[] foox = new double[1000000];
		double[] fooy = new double[1000000];
		
		// gather data
		rin = new ReadInUtil(params.sget("Input File (tab del.)"));
		rin.getData(new int[]{params.iget("Time Col. Num."),params.iget("Metric Col. Num.")},1,foox,fooy);
		
		while(foox[count] > 0){
			count++;
		}
		
		// copy data without zero padding
		xvals = CopyUtil.copyArray(foox,count);		
		yvals = CopyUtil.invertAndScaleArray(fooy,count,1);
		
		// fit data
//	//	fitter = new FitUtil(count);
//		//double temp[] = fitter.fit(xvals, yvals, 1,false);
//		
//		// fit params
//		//m    = temp[0];
//		b    = temp[1];
//		chi2 = temp[4];
//		
//		System.out.println(m);
//		System.out.println(b);
//		System.out.println(chi2);
		
		Job.animate();
		
		return;
	}
	
	

}

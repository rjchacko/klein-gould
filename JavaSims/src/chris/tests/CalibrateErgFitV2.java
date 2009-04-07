package chris.tests;

import java.io.File;

import scikit.jobs.Control;
import scikit.jobs.Job;
import scikit.jobs.Simulation;
import scikit.jobs.params.DirectoryValue;
import chris.util.CopyUtil;
import chris.util.FitUtil;
import chris.util.PrintUtil;
import chris.util.ReadInUtil;

public class CalibrateErgFitV2 extends Simulation{

	private ReadInUtil rin;
	private FitUtil fitter;
	private static String[] fins = new String[]{"0p","1p","10p","15p","16p","17p","18p","20p","30p","50p"};
	//private static String[] fins = new String[]{"16p"};
	private String din, fin, fout;
	private double foox[], fooy[];
	
	
	public static void main(String[] args) {
		new Control(new CalibrateErgFitV2(), "Calibrate Ergodic Fitter");
	}
	
	public void animate() {

		return;
	}

	
	public void clear() {

		return;
	}

	
	public void load(Control c) {
		params.add("Input File (tab del.)",new DirectoryValue("/Users/cserino/Research/Data2/ErgPhaseDiagram/NearestNeighbor/"));
		params.add("Time Col. Num.", (int) 1);
		params.add("Metric Col. Num.", (int) 8);
		params.add("File Number");
		params.set("File Number",fins[0]);
	}

	public void run() {
		
		din  = params.sget("Input File (tab del.)");
		fout = din + File.separator + "FitData.txt";
		
		PrintUtil.printlnToFile(fout,"File","slope","intercept","chi^2");
		
		
		for(int jj = 0 ; jj < fins.length ; jj++){
		
			params.set("File Number",fins[jj]);
			Job.animate();	
			
			
			int count = 0;

			fin = din+ File.separator + "StressData_" + fins[jj] + ".txt";
			
			foox = new double[1000000];
			fooy = new double[1000000];
			
			// gather data
			rin = new ReadInUtil(fin);
			rin.getData(new int[]{params.iget("Time Col. Num."),params.iget("Metric Col. Num.")},1,foox,fooy);

			if(foox[foox.length - 1] == 0){
				while(foox[count] > 0){
					count++;
				}
				// copy data without zero padding
				foox = CopyUtil.copyArray(foox,count);		
				fooy = CopyUtil.invertAndScaleArray(fooy,count,1);
			}

			// fit data
			fitter = new FitUtil(count);
			double temp[] = fitter.fit(foox, fooy, 1, false);

			// fit params
			PrintUtil.printlnToFile(fout,fins[jj],temp[0],temp[1],temp[4]);
	
		}
		
		params.set("File Number","Done");;
		Job.animate();
		
		return;
	}
	
	

}

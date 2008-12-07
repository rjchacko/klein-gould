package chris.RandomWalker;

import java.io.File;
import java.text.DecimalFormat;
import java.util.Random;

import chris.util.PrintUtil;

import scikit.jobs.Control;
import scikit.jobs.Job;
import scikit.jobs.Simulation;
import scikit.jobs.params.DirectoryValue;
import scikit.jobs.params.DoubleValue;
import scikit.jobs.params.StringValue;


public class GeometricStepsLoop extends Simulation{

	private int seed, walkers, steps, pdf[], bins;
	private double lambda, drmax, binsize;
	private String outdir, fout, pfile;
	private Random rand;
	private boolean r2o;
	
	DecimalFormat fmt = new DecimalFormat("0.00");
	
	private static final double twoPI = 8*Math.atan(1);
	
	public static void main(String[] args) {
		new Control(new GeometricStepsLoop(), "Parameters");
	}
	
	public void load(Control c) {
	
		params.add("Data Directory",new DirectoryValue("/Users/cserino/CurrentSemester/PY542/Final"));
		params.add("File Name", new StringValue("foo"));
		params.add("Random Seed",(int) 0);
		params.add("Bins", (int) 500);
		params.add("Walkers", (int) 10);
		params.add("Steps", (int) 10);
		params.add("lambda", new DoubleValue(0.01, 0., 1.));
		params.add("Status");
	}

	public void run() {
	
		// Step up parameters for the simulation
		
		while(lambda < 1){
		
			setup();

			params.set("Status", "Running");
			// Perform random walks
			for (int jj = 0 ;jj < walkers ; jj++){
				updatePDF(walk());
				if(jj%10000 == 0){
					params.set("Status", jj);
					Job.animate();
				}
			}

			savedata();
			params.set("Status", "Done");
			Job.animate();
		    
			params.set("Random Seed", params.iget("Random Seed")+1);
			params.set("lambda", params.fget("lambda")+0.01);
			
		}
	
	}
	
	private void setup(){
		
		String tmp;
		
		seed    = params.iget("Random Seed");
		bins    = params.iget("Bins");
		walkers = params.iget("Walkers");
		steps   = params.iget("Steps");
		lambda  = params.fget("lambda");
		outdir  = params.sget("Data Directory");
		tmp     = params.sget("File Name");
		
//		local = new double[walkers][2];
		fout  = outdir + File.separator + tmp + "_" + fmt.format(lambda) +".txt";
		pfile = outdir + File.separator + "Params_" + tmp + "_" + fmt.format(lambda)  + ".txt";
		rand  = new Random(seed); 
		pdf   = new int[bins];
		
		PrintUtil.printlnToFile(pfile, params.toString());
		
		drmax = 1;
		for(int jj = 0 ; jj < bins ; jj++){	// do we ever have fewer bins than steps????
			pdf[jj] = 0;
			if (jj < steps){
				drmax = drmax*lambda;
			}
		}
		
		
		if(lambda != 1){
			drmax   = (lambda - drmax)/(1-lambda);
			if( drmax < 1){
				binsize = 2*drmax/bins; 
				r2o = false;
				PrintUtil.printlnToFile(pfile, "bin 0 corresponds to r = 1 - dr_max = ", 1-drmax);
				PrintUtil.printlnToFile(pfile, "binsize = ", binsize);
			}
			else{
				binsize = (1.+drmax)/bins; 
				r2o = true;
				PrintUtil.printlnToFile(pfile, "bin 0 corresponds to r = 0 = ", 0.);
				PrintUtil.printlnToFile(pfile, "binsize = ", binsize);
			}
		}
		else{
			drmax = steps;
			binsize = drmax/bins;
			r2o = true;
			PrintUtil.printlnToFile(pfile, "bin 0 corresponds to r = 0 = ", 0.);
			PrintUtil.printlnToFile(pfile, "binsize = ", binsize);
		}

		return;
	}
	
	private double[] walk(){
		
		double stepsize = 1;
		double[] rr = new double[2];
		double theta, costheta, sintheta;
		
		rr[0] = 0;
		rr[1] = 0;
		
		for (int jj = 0 ; jj < steps ; jj++){
			/*
			 * This is probably slow!
			 */
			theta = nextAngle();
			costheta = Math.cos(theta);
			sintheta = Math.sin(theta);
			rr[0] += stepsize*costheta;
			rr[1] += stepsize*sintheta;
			
			stepsize = stepsize*lambda;
		}
		
		return rr;
	}
	
	private double nextAngle(){
		/*
		 * might be faster to just generate cosines (or sines)
		 */
		return twoPI*rand.nextDouble();
	}
	
	private void updatePDF(double[] rr){
		
		/*
		 * SHIFT IF drmax < 1
		 * 
		 */
		if(r2o){
			pdf[(int)(Math.sqrt(rr[0]*rr[0] + rr[1]*rr[1])/binsize)]++;	
		}
		else{
			pdf[(int)((Math.sqrt(rr[0]*rr[0] + rr[1]*rr[1])-1.+drmax)/binsize)]++;	
		}
		return;
	}
	
	private void savedata(){

		if(r2o){
			PrintUtil.printWalkerData(fout,pdf,bins,0,binsize);
		}
		else{
			PrintUtil.printWalkerData(fout,pdf,bins,1-drmax,binsize);
		}
		return;
	}
	
	public void animate() {
	
		return;
	}

	public void clear() {
	
		return;
	}

	


}

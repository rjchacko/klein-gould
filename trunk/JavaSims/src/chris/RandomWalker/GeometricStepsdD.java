package chris.RandomWalker;

import java.io.File;
import java.util.Random;

import chris.util.PrintUtil;

import scikit.jobs.Control;
import scikit.jobs.Job;
import scikit.jobs.Simulation;
import scikit.jobs.params.DirectoryValue;
import scikit.jobs.params.DoubleValue;
import scikit.jobs.params.StringValue;


public class GeometricStepsdD extends Simulation{

	private int seed, walkers, steps, pdf[], bins, dim;
	private double lambda, drmax, binsize;
	//private double local[][];
	private String outdir, fout, pfile;
	private Random rand;
	private boolean r2o;
		
	public static void main(String[] args) {
		new Control(new GeometricStepsdD(), "Parameters");
	}
	
	public void load(Control c) {
	
		params.add("Data Directory",new DirectoryValue("/Users/cserino/CurrentSemester/PY542/Final"));
		params.add("File Name", new StringValue("foo"));
		params.add("Random Seed",(int) 0);
		params.add("Dimension", (int) 2);
		params.add("Bins", (int) 500);
		params.add("Walkers", (int) 100000000);
		params.add("Steps", (int) 20);
		params.add("lambda", new DoubleValue(0.5, 0., 1.));
		params.add("Status");
	}

	public void run() {
	
		// Step up parameters for the simulation
		setup();

		params.set("Status", "Running");
		// Perform random walks
		for (int jj = 0 ;jj < walkers ; jj++){
//			local[jj] = walk();
//			updatePDF(local[jj]);
			updatePDF(walk());
			if(jj%10000 == 0){
				params.set("Status", jj);
				Job.animate();
			}
		}

		savedata();
		params.set("Status", "Done");
		Job.animate();
		
		
	
	}
	
	private void setup(){
		
		String tmp;
		
		seed    = params.iget("Random Seed");
		dim     = params.iget("Dimension");
		bins    = params.iget("Bins");
		walkers = params.iget("Walkers");
		steps   = params.iget("Steps");
		lambda  = params.fget("lambda");
		outdir  = params.sget("Data Directory");
		tmp     = params.sget("File Name");
		
//		local = new double[walkers][2];
		fout  = outdir + File.separator + tmp + ".txt";
		pfile = outdir + File.separator + "Params_" + tmp + ".txt";
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
		double[] rr = new double[dim];
		
		for (int jj = 0 ; jj < dim ; jj++){
			rr[jj] = 0;
		}

		
		for (int jj = 0 ; jj < steps ; jj++){

			double[] tmp = nextStep(stepsize);

			for (int kk = 0 ; kk < dim ; kk++){
				rr[kk] += tmp[kk];
			}
			
			stepsize = stepsize*lambda;
		}
		
		return rr;
	}
	
	private double[] nextStep(double norm){
		double[] ret = new double[dim];
		double tmpnorm = 0;
		
		for (int jj = 0 ; jj < dim ; jj++){
			ret[jj] = rand.nextDouble()-0.5;
			tmpnorm += ret[jj]*ret[jj];
		}
		tmpnorm = Math.sqrt(tmpnorm);
		for (int jj = 0 ; jj < dim ; jj++){
			ret[jj] = norm*ret[jj]/tmpnorm;
		}
		
		return ret;
	}
	
	private void updatePDF(double[] rr){
		
		/*
		 * SHIFT IF drmax < 1
		 * 
		 */
		
		double rs = 0;
		for (int jj = 0 ; jj < dim ; jj++){
			rs += rr[jj]*rr[jj];
		}
		rs = Math.sqrt(rs);
		
		if(r2o){
			pdf[(int)(rs/binsize)]++;	
		}
		else{
			pdf[(int)(rs/binsize)]++;	
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

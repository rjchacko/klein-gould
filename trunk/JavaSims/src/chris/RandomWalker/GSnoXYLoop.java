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


public class GSnoXYLoop extends Simulation{

	private int seed, walkers, steps, pdf[][], bins;
	private double lambda, binsize;
	private String outdir, fout, pfile;
	private Random rand;
	
	private static final double twoPI = 8*Math.atan(1);
	
	public static void main(String[] args) {
		new Control(new GSnoXYLoop(), "Parameters");
	}
	
	public void load(Control c) {
	
		params.add("Data Directory",new DirectoryValue("/Users/cserino/CurrentSemester/PY542/Final/Data"));
		params.add("File Name", new StringValue("foo"));
		params.add("Random Seed",(int) 0);
		params.add("Bins", (int) 1000);
		params.add("Walkers", (int) 1E7);
		params.add("Steps", (int) 50);
		params.add("lambda", new DoubleValue(0.75, 0.5, 1.));
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
		pdf   = new int[bins][bins];
		
		PrintUtil.printlnToFile(pfile, params.toString());
		
		for(int jj = 0 ; jj < bins ; jj++){	
			for(int kk = 0 ; kk < bins ; kk++){
				pdf[jj][kk] = 0;
			}
		}
		
		PrintUtil.printlnToFile(pfile, "bin 0 corresponds to r = 0 = ", 0.);
		PrintUtil.printlnToFile(pfile, "binsize = ", binsize);

		return;
	}
	
	private double[] walk(){
		
		double stepsize = 1;
		double[] rr = new double[2];
		double theta, costheta, sintheta;
		
		rr[0] = 0;
		rr[1] = 0;
		
		for (int jj = 0 ; jj < steps ; jj++){
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
		
		// have *bins* number of bins (for x and y separately) which we want to position 
		// evenly about x=y=0. 
		
		/* (-0.5,0.5)	(0.5,0.5)
		 * 		xxxxxxxxxxxxx
		 * 		x			x	
		 * 		x			x	
		 * 		x			x
		 * 		xxxxxxxxxxxxx		
		 * (-0.5,-0.5)	(0.5,-0.5)	
		 * 
		 */		
		
		if(Math.abs(rr[0])<0.5 && Math.abs(rr[1])<0.5 ) pdf[(int)((rr[0]+0.5)*bins)][(int)((rr[1]+0.5)*bins)]++;
		
		return;
	}
	
	private void savedata(){
		
		PrintUtil.printlnToFile(fout,"x","y","pdf");
		PrintUtil.printWalkerData(fout,pdf,bins);

		return;
	}
	
	public void animate() {
	
		return;
	}

	public void clear() {
	
		return;
	}

	


}

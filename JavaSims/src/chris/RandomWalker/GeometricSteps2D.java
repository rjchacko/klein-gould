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


public class GeometricSteps2D extends Simulation{

	private int seed, walkers, steps, pdfR[], pdfX[], pdfY[], bins;
	private double lambda, drmax, binsize, gsum;
	//private double local[][];
	private String outdir, foutR, foutX, foutY, pfile;
	private Random rand;
	private boolean r2o;
	
	private static final double twoPI = 8*Math.atan(1);
	
	public static void main(String[] args) {
		new Control(new GeometricSteps2D(), "Parameters");
	}
	
	public void load(Control c) {
	
		params.add("Data Directory",new DirectoryValue("/Users/cserino/Desktop"));
		params.add("File Name", new StringValue("default"));
		params.add("Random Seed",(int) 0);
		params.add("Resolution", (double) 0.01);
		params.add("Walkers", (int) 1e8);
		params.add("Steps", (int) 10);
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
		binsize = params.fget("Resolution");
		walkers = params.iget("Walkers");
		steps   = params.iget("Steps");
		lambda   = params.fget("lambda");
		outdir  = params.sget("Data Directory");
		tmp     = params.sget("File Name");
		pfile   = outdir + File.separator + "Params_" + tmp + ".txt";
		
		PrintUtil.printlnToFile(pfile, params.toString());
		
		drmax = 1;
		for(int jj = 0 ; jj < bins ; jj++){	
			pdfR[jj] = 0;
		}
		
		if(lambda != 1){
			gsum = 1./(1.-lambda);
			if(lambda < 0.5){
				r2o = false;
				bins = (int)((2*lambda/(1-lambda))/binsize);  
			}
			else{
				r2o  = true;
				bins = (int)((1+(lambda/(1-lambda)))/binsize);  
			}
		}
		else{
			r2o = true;
			bins = (int)(steps/binsize);
			gsum = steps;
		}
		PrintUtil.printlnToFile(pfile, "binsize = ", binsize);

		foutR = outdir + File.separator + tmp + "_R.txt";
		foutX = outdir + File.separator + tmp + "_X.txt";
		foutY = outdir + File.separator + tmp + "_Y.txt";
		rand  = new Random(seed); 
		pdfR  = new int[bins];
		pdfX  = new int[2*bins];
		pdfY  = new int[2*bins];
		
		
		
		
		
		
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
			pdfR[(int)(Math.sqrt(rr[0]*rr[0] + rr[1]*rr[1])/binsize)]++;	
			pdfX[(int)((rr[0]+gsum)/binsize)]++;
			pdfY[(int)((rr[1]+gsum)/binsize)]++;
		}
		else{
			pdfR[(int)((Math.sqrt(rr[0]*rr[0] + rr[1]*rr[1])-1.+drmax)/binsize)]++;	
		}
		return;
	}
	
	private void savedata(){

		if(r2o){
			PrintUtil.printWalkerData(foutR,pdfR,bins,0,binsize);
			PrintUtil.printWalkerData(foutX,pdfX,2*bins,-gsum,binsize);
			PrintUtil.printWalkerData(foutY,pdfY,2*bins,-gsum,binsize);
		}
		else{
			PrintUtil.printWalkerData(foutR,pdfR,bins,1-drmax,binsize);
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

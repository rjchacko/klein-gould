package chris.RandomWalker;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.text.DecimalFormat;
import java.util.Random;

import chris.util.PrintUtil;

import scikit.jobs.Control;
import scikit.jobs.Job;
import scikit.jobs.Simulation;
import scikit.jobs.params.DirectoryValue;
import scikit.jobs.params.DoubleValue;
import scikit.jobs.params.StringValue;


public class GSconditional2D extends Simulation{

	private int seed, walkers, steps;
	private double lambda, gsum;
	//private double local[][];
	private String outdir, pfile;
	private Random rand;
	private DecimalFormat fmt = new DecimalFormat("0.00");
	double[][] pdfR = new double[50][50];
	double[][] pdfX = new double[50][50];
	private static final double twoPI = 8*Math.atan(1);
	
	public static void main(String[] args) {
		new Control(new GSconditional2D(), "Parameters");
	}
	
	public void load(Control c) {
	
		params.add("Data Directory",new DirectoryValue("/Users/cserino/Desktop/"));
		params.add("File Name", new StringValue("foo"));
		params.add("Random Seed",(int) 0);
		params.add("Walkers", (int) 1e6);
		params.add("Steps", (int) 10);
		params.add("lambda", new DoubleValue(0.5, 0.5, 1.));
		params.add("Status");
	}

	public void run() {
	
		// Step up parameters for the simulation
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
		params.set("Status", "Saving");
		Job.animate();
		savedata();
		
		params.set("Status", "Done");
		Job.animate();
		
		
	
	}
	
	private void setup(){
		
		String tmp;
		
		seed    = params.iget("Random Seed");
		walkers = params.iget("Walkers");
		steps   = params.iget("Steps");
		lambda  = params.fget("lambda");
		outdir  = params.sget("Data Directory");
		tmp     = params.sget("File Name");
		
		pfile = outdir + File.separator + "Params_" + tmp + ".txt";
		rand  = new Random(seed); 

		PrintUtil.printlnToFile(pfile, params.toString());

		gsum = 1 / ( 1 - lambda );
		
		for(int jj = 0 ; jj < 50 ; jj++){
			for(int kk = 0 ; kk < 50 ; kk++){
				pdfR[jj][kk] = 0;
				pdfX[jj][kk] = 0;
			}
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

		double thisR;

		thisR = Math.sqrt(sqr(rr[0])+ sqr(rr[1]));
		// pdfR [ R ] [ x ]
		pdfR[(int)(50*thisR/gsum)][(int)(25*(rr[0]+gsum)/gsum)]++;
		// pdfX [ x ] [ R ]
		pdfX[(int)(25*(rr[0]+gsum)/gsum)][(int)(50*thisR/gsum)]++;
		                                    
		return;
	}
	
	private void savedata(){

		String xname, rname;
		
		for(int jj = 0 ; jj < 50 ; jj++){

			if(((jj/25)-1) < 0){
				xname = outdir + File.separator + "X_neg"+ fmt.format(-((jj/25.)-1.)*gsum)+".txt";
			}
			else{
				xname = outdir + File.separator + "X_"+ fmt.format(((jj/25.)-1.)*gsum)+".txt";

			}			
			rname = outdir + File.separator + "R_"+ fmt.format(jj*gsum/50)+".txt";

			try{
				File fx = new File(xname);
				File fr = new File(rname);
				PrintWriter pwx = new PrintWriter(new FileWriter(fx, true), true);
				PrintWriter pwr = new PrintWriter(new FileWriter(fr, true), true);
				pwx.print("r");
				pwx.print("\t");
				pwx.print("freq");
				pwx.println();
				pwr.print("x");
				pwr.print("\t");
				pwr.print("freq");
				pwr.println();
				for(int kk = 0 ; kk < 50 ; kk++){
					pwx.print(kk*gsum/50);
					pwx.print("\t");
					pwx.print(pdfX[jj][kk]);
					pwx.println();
					pwr.print(((kk/25.)-1.)*gsum);
					pwr.print("\t");
					pwr.print(pdfR[jj][kk]);
					pwr.println();
				}
				pwx.close();
				pwr.close();
			}
			catch (IOException ex){
				ex.printStackTrace();
			}
		}

		return;
	}
	
	public void animate() {
	
		return;
	}

	public void clear() {
	
		return;
	}
	
	private double sqr(double x){
		
		return x*x;
	}

}

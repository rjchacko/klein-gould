package chris.RandomWalker;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.text.DecimalFormat;
import java.util.Random;

import scikit.jobs.Control;
import scikit.jobs.Job;
import scikit.jobs.Simulation;
import scikit.jobs.params.DirectoryValue;
import scikit.jobs.params.StringValue;
import chris.util.PrintUtil;


public class PowerLawSteps2DLoop extends Simulation{

	private int seed, walkers, steps, pdf[];
	private double alpha, binsize, alphaMx, da, slength[];
	private String outdir, fout, pfile, bname;
	private Random rand;
	
	DecimalFormat fmt = new DecimalFormat("0.00");
	
	private static final double twoPI = 8*Math.atan(1);
	
	public static void main(String[] args) {
		new Control(new PowerLawSteps2DLoop(), "Parameters");
	}
	
	public void load(Control c) {
	
		params.add("Data Directory",new DirectoryValue("/Users/cserino/Desktop"));
		params.add("File Name", new StringValue("foo"));
		params.add("Random Seed",(int) 0);
		params.add("Bin Size", (double) 0.1);
		params.add("Walkers", (int) 10000000);
		params.add("alpha Min", (double) 0);	// N**(1-2a)
		params.add("alpha Max", (double) 1.5);	// Sum converges to \zeta(2a-1)
		params.add("d(alpha)", (double) 0.1);
		params.add("Status");
	}

	public void run() {
	
		// Step up parameters for the simulation
		alpha   = params.fget("alpha Min");
		alphaMx = params.fget("alpha Max");
		da      = params.fget("d(alpha)");
		binsize = params.fget("Bin Size");
		walkers = params.iget("Walkers");
		outdir  = params.sget("Data Directory");
		bname   = params.sget("File Name");
		
		while(alpha <= alphaMx){
		
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
			alpha += da;
		}
	
	}
	
	private void setup(){

		int Nbins;
		double zeta;
		
		seed    = params.iget("Random Seed");		
		fout    = outdir + File.separator + bname + "_" + fmt.format(alpha) +".txt";
		pfile   = outdir + File.separator + "Params_" + bname + "_" + fmt.format(alpha)  + ".txt";
		rand    = new Random(seed); 
		steps = (int)(Math.ceil(Math.pow(binsize,(1./(1.-2*alpha))))+5);
		if(steps > 100) steps = 100;
		slength = new double[steps];
		
		zeta = 0;
		for(int jj = 1 ; jj <= steps ; jj++ ){
			slength[jj-1] = Math.pow(jj,1-2*alpha);
			zeta += slength[jj-1];
		}
		Nbins = (int)(Math.ceil(zeta/binsize)+1);
		pdf   = new int[Nbins];
		
		PrintUtil.printlnToFile(pfile, params.toString());
		PrintUtil.printlnToFile(pfile, steps);
		return;
	}
	
	private double[] walk(){

		double[] rr = new double[2];
		double theta, costheta, sintheta;
		
		rr[0] = 0;
		rr[1] = 0;
		
		for (int jj = 0 ; jj < steps ; jj++){
			theta = nextAngle();
			costheta = Math.cos(theta);
			sintheta = Math.sin(theta);
			rr[0] += slength[jj]*costheta;
			rr[1] += slength[jj]*sintheta;
		}
		
		return rr;
	}
	
	private double nextAngle(){

		return twoPI*rand.nextDouble();
	}
	
	private void updatePDF(double[] rr){
		
		double r = Math.sqrt(rr[0]*rr[0]+rr[1]*rr[1]);
		pdf[(int)(r/binsize)]++;
		
		return;
	}
	
	private void savedata(){
		
		try{
			File file = new File(fout);
			PrintWriter pw = new PrintWriter(new FileWriter(file, true), true);
			for (int ii = 0 ; ii < pdf.length ; ii++){				
				pw.print(ii);
				pw.print("\t");
				pw.print(pdf[ii]);
				pw.println();
			}
			pw.close();
		}
		catch (IOException ex){
			ex.printStackTrace();
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

package rachele.ising.dim2.apps;

import java.awt.Color;
import rachele.ising.dim2.*;
import rachele.util.FourierTransformer;
import scikit.dataset.Accumulator;
import scikit.dataset.Function;
//import scikit.dataset.PointSet;
import scikit.graphics.dim2.Grid;
import scikit.graphics.dim2.Plot;
import scikit.jobs.Control;
import scikit.jobs.Job;
import scikit.jobs.Simulation;
import scikit.jobs.params.ChoiceValue;
import scikit.jobs.params.DoubleValue;
import static java.lang.Math.exp;
//import static java.lang.Math.floor;


/**
* 
* Langevin field theory application.  Plots average of Svt and varience of the
* order parameter for max time given by the user.
* Also plots predicted value of linear theory.
* This should work now for critical and off-critical
* quenches.  For off-critical quenches, use with 1D App 
* to find out the proper input magnetization.  Start it with no noise and it
* will quickly find the uniform background configuration.  Input this
* into this program for magnetization.  Use Random Guassian init conditions.
* 
*/
public class LinearTheoryApp extends Simulation{
    Grid grid = new Grid("Phi(x)");
	Plot structurePeak = new Plot("S v t");
	//Plot structurePeakH = new Plot("Hor Structure factor");
	Plot sfPeakBoth = new Plot("Both Structure factors");
	Plot variance = new Plot("Variance");
	Plot sfVsKR = new Plot("Structure Function");
	Plot sfSlice = new Plot("SF Slice");
	//StructureFactor sf;
    FourierTransformer ft;
	IsingField2D ising;
    int kR1int,kR2int;
	Accumulator peak1AveAcc;
	Accumulator peak1Acc;
    Accumulator varianceAcc;
    Accumulator meanPhiAcc;
    double phi0;
    boolean langevin_dynamics;
    
	public static void main(String[] args) {
		new Control(new LinearTheoryApp(), "Ising Field");
	}

	public void load(Control c) {
		c.frameTogether("Plots", grid, sfVsKR, structurePeak, variance);
		//c.frame(sfSlice);
		params.addm("Zoom", new ChoiceValue("Yes", "No"));
		params.addm("Interaction", new ChoiceValue("Square", "Circle"));
		params.addm("Dynamics?", new ChoiceValue("Langevin Conserve M", "Langevin No M Convervation"));
		params.add("Init Conditions", new ChoiceValue("Random Gaussian", 
				"Artificial Stripe 3", "Read From File", "Constant" ));
		params.addm("Approx", new ChoiceValue("Exact Stable",
				"Avoid Boundaries", "Exact SemiStable", "Exact", "Linear",  "Phi4"));
		params.add("Dynamics", new ChoiceValue( "Glauber", "Langevin"));
		params.addm("Noise", new DoubleValue(1.0, 0.0, 1.0).withSlider());
		params.addm("Horizontal Slice", new DoubleValue(0.5, 0, 0.9999).withSlider());
		params.addm("Vertical Slice", new DoubleValue(0.5, 0, 0.9999).withSlider());
		params.addm("T", 0.04);
		params.addm("H", 0.8);
		params.add("Magnetization", 0.760138297217);
		params.addm("dt", 0.001);
		params.addm("J", -1.0);
		params.addm("R", 200000.0);
		params.add("L/R", 2.782608696);
		params.add("R/dx", 50.0);
		params.add("kR bin-width", 0.1);
		params.add("Random seed", 0);
		params.add("Max Time", 3.0);
		params.add("Time");
		params.add("Reps");
		params.add("Mean Phi");
		params.add("Lp");
		
	}
	
	public void animate() {
		ising.readParams(params);
		params.set("Time", ising.time());
		if (params.sget("Zoom").equals("Yes")) 	grid.setAutoScale();
		else grid.setScale(-1, 1);
		structurePeak.setAutoScale(true);
		structurePeak.setLogScale(false, true);
		grid.registerData(ising.Lp, ising.Lp, ising.phi);

			structurePeak.registerLines("Structure theory", new Function(0, params.fget("Max Time")) {
				public double eval(double t) {
					double kR = 2.0*Math.PI*kR1int*ising.R/ising.L;
					double pot = (kR == 0) ? 1 : Math.sin(kR)/kR; 
					double D = (-pot-ising.T/(1.0-phi0*phi0))*(1-phi0*phi0)/ising.T;
					return  (exp(2*t*D)*(1 + 1/D)-1/D);
				}
			}, Color.BLUE);
			structurePeak.registerLines("Structure 2 theory", new Function(0, params.fget("Max Time")) {
				public double eval(double t) {
					double kR = 2.0*Math.PI*kR1int*ising.R/ising.L;
					double pot = (kR == 0) ? 1 : Math.sin(kR)/kR; 
					double D = (-pot-ising.T/(1.0-phi0*phi0))*(1-phi0*phi0)/ising.T;
					return  (exp(2*t*D));
				}
			}, Color.GREEN);
		structurePeak.registerLines("First Peak", peak1AveAcc, Color.BLACK);
		structurePeak.registerLines("First Peak Instance", peak1Acc, Color.RED);
		variance.registerLines("Variance", varianceAcc, Color.BLACK);
	}

	public void clear() {
	}

	public void run() {
		ising = new IsingField2D(params);
		ft = new FourierTransformer(ising.Lp);
		peak1AveAcc = new Accumulator(ising.dt*2);
		peak1AveAcc.enableErrorBars(true);
		peak1Acc = new Accumulator(ising.dt*2);
	    varianceAcc = new Accumulator(ising.dt);
	    meanPhiAcc = new Accumulator(ising.dt);
		varianceAcc.clear();
		meanPhiAcc.clear();
	    phi0 = params.fget("Magnetization");
		kR1int = ising.findKRSquareInt(IsingField2D.KRsquare);
		kR2int = ising.findKRSquareInt(IsingField2D.KRsquare2);
		double maxTime = params.fget("Max Time");
		if(params.sget("Dynamics") == "Langevin") langevin_dynamics = true;
		else langevin_dynamics = false;

		int reps = 0;
//		double abortTimeSum = 0;
		
		while (true) {
			peak1Acc.clear();
			ising.randomizeField(phi0);
			ising.restartClock();
			boolean abort = false;
			collect();
			Job.animate();
			params.set("dt", 0.01);
			while(ising.time()< maxTime & abort == false){
				//if(ising.time()>4.5) params.set("dt", 0.001);
				params.set("Mean Phi", ising.mean(ising.phi));
				if(langevin_dynamics) ising.simulateSimple();
				else ising.simulateGlauber();
//				abort = ising.simulateUnstable();	
//				if (abort){
//					System.out.println("abort at " + ising.time());
//					abortTimeSum += ising.time();
//				}else collect();
				collect();
				Job.animate();
			}
		reps += 1;
		params.set("Reps", reps);
//		System.out.println(abort);
//		double aveAbortTime = abortTimeSum/(double)reps;
//		System.out.println("ave abort time = " + aveAbortTime);
		}
	}
	
	public void collect(){
		double [] sf = ft.find2DSF(ising.phi, ising.L);
		peak1AveAcc.accum(ising.time(),sf[kR1int]);
		peak1Acc.accum(ising.time(),sf[kR1int]);
		varianceAcc.accum(ising.time(), ising.phiVariance());
		meanPhiAcc.accum(ising.time(),ising.mean(ising.phi));		
	}
	
	public double findMeanPhi(){
		for (double t = 0.0; t < 50.0; t = t + params.fget("dt")){
			ising.simulate();
			System.out.println("init for density: time = " + ising.time());
		}
		return ising.mean(ising.phi);
	}
}

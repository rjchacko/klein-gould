package rachele.ising.dim2.apps;

import java.awt.Color;
import rachele.ising.dim2.*;
import rachele.util.FileUtil;
import rachele.util.FourierTransformer;
import scikit.dataset.Accumulator;
import scikit.dataset.Function;
import scikit.graphics.dim2.Grid;
import scikit.graphics.dim2.Plot;
import scikit.jobs.Control;
import scikit.jobs.Job;
import scikit.jobs.Simulation;
import scikit.jobs.params.ChoiceValue;
import scikit.jobs.params.DoubleValue;
import static java.lang.Math.*;

/**
*Creates a plot of Structure factor vs kR for a 
*critical quench from disorder to below the critical temperature.
*Writes to file with error bars.
*/
public class LinearTheorySkApp extends Simulation{
    Grid grid = new Grid("Phi(x)");
	Plot sfSlice = new Plot("SF Slice");
	Accumulator sfAcc,sfInstance;
	FourierTransformer ft;
    IsingField2D ising;
    
    
	public static void main(String[] args) {
		new Control(new LinearTheorySkApp(), "Ising Field");
	}

	public void load(Control c) {
		c.frame(grid);
		c.frame(sfSlice);
		params.addm("Interaction", new ChoiceValue("Square", "Circle"));
		params.addm("Dynamics?", new ChoiceValue("Langevin Conserve M", "Langevin No M Convervation"));
		params.add("Init Conditions", new ChoiceValue("Random Gaussian", 
				"Artificial Stripe 3", "Read From File", "Constant" ));
		params.addm("Approx", new ChoiceValue("Exact Stable",
				"Avoid Boundaries", "Exact SemiStable", "Exact", "Linear",  "Phi4"));
		params.addm("Noise", new DoubleValue(1.0, 0.0, 1.0).withSlider());
		params.addm("Horizontal Slice", new DoubleValue(0.5, 0, 0.9999).withSlider());
		params.addm("Vertical Slice", new DoubleValue(0.5, 0, 0.9999).withSlider());
		params.addm("T", 0.04);
		params.addm("H", 0.0);
		params.add("Magnetization", 0.0);
		params.addm("dt", 0.001);
		params.addm("J", -1.0);
		params.addm("R", 1000000.0);
		params.add("L/R", 20.0);
		params.add("R/dx", 4.0);
		params.add("kR bin-width", 0.1);
		params.add("Random seed", 0);
		params.add("Max Time", 0.125);
		params.add("Time");
		params.add("Reps");
		params.add("Mean Phi");
		params.add("Lp");
		
	}
	
	public void animate() {
		ising.readParams(params);
		params.set("Time", ising.time());
		grid.registerData(ising.Lp, ising.Lp, ising.phi);
		sfSlice.registerLines("slice2", sfAcc, Color.BLUE);
		sfSlice.registerLines("slice3", sfInstance, Color.YELLOW);
		double maxValue = PI*ising.Lp*ising.R/ising.L;
		sfSlice.registerLines("Structure theory", new Function(0, maxValue) {
			public double eval(double kR) {
				double pot = (kR == 0) ? 1 : Math.sin(kR)/kR; 
				double D = -pot/ising.T-1;
				return  exp(2*ising.time()*D);//*(1 + 1/D)-1/D;	
			}
		}, Color.RED);
		sfSlice.registerLines("Structure theory noise", new Function(0, maxValue) {
			public double eval(double kR) {
				double pot = (kR == 0) ? 1 : Math.sin(kR)/kR; 
				double D = -pot/ising.T-1;
				return  (exp(2*ising.time()*D)*(1 + 1/D)-1/D);	
			}
		}, Color.GREEN);
	}

	public void clear() {
	}

	public void run() {
		sfAcc = new Accumulator(); sfAcc.enableErrorBars(true);
		sfInstance = new Accumulator();
		ising = new IsingField2D(params);
		ft = new FourierTransformer(ising.Lp);
		//double binWidth = 2.0*PI*ising.R/ising.Lp;
		double maxTime = params.fget("Max Time");
		int reps = 0;
		while (true) {
			sfInstance.clear();
			ising.randomizeField(0);
			ising.restartClock();
			while(ising.time() <= maxTime){
				ising.simulateSimple();
				Job.animate();
			}
			//have to input dx into FFT. I don't know why, but its true
			double [] sf2 = ft.find2DSF(ising.phi, ising.L);
			for (int i = 0; i < ising.Lp/2; i++){
				double kRValue = 2*PI*i*ising.R/ising.L;
				sfAcc.accum(kRValue, sf2[i]);
				sfInstance.accum(kRValue, sf2[i]);
			}

			reps += 1;
			params.set("Reps", reps);
			if(reps%2 == 0) writeToFile(ising.dt, reps);
		}
	}
	
	private void writeToFile(double recordTime1, int reps){
		String message1 = "#Field theory data: S vs k for one or more times. Disorder to Order Early times.";
		StringBuffer sb = new StringBuffer();
		sb.append("# record time = ");
		sb.append(recordTime1);
		sb.append(", rep no = ");
		sb.append(reps);
		String message2 = sb.toString();
		String fileName = "../../../research/javaData/stripeToClumpInvestigation/monteCarloData/squareResults/svkCompare4/f";
		FileUtil.initFile(fileName, params, message1, message2);
		FileUtil.printAccumToFile(fileName, sfAcc);
		System.out.println("file written for rep " + reps);
	}

}

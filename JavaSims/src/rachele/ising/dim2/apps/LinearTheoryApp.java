package rachele.ising.dim2.apps;

import java.awt.Color;

import rachele.ising.dim2.*;
import scikit.dataset.Accumulator;
import scikit.graphics.dim2.Grid;
import scikit.graphics.dim2.Plot;
import scikit.jobs.Control;
import scikit.jobs.Job;
import scikit.jobs.Simulation;
import scikit.jobs.params.ChoiceValue;
import scikit.jobs.params.DoubleValue;
import static java.lang.Math.PI;
import static java.lang.Math.abs;
import static java.lang.Math.floor;
import static scikit.numerics.Math2.j1;
import static scikit.util.Utilities.frameTogether;
import static java.lang.Math.*;

public class LinearTheoryApp extends Simulation{
    Grid grid = new Grid("Phi(x)");
	Plot structurePeakV = new Plot("Ver Structure Factor");
	Plot structurePeakH = new Plot("Hor Structure factor");
	Plot sfPeakBoth = new Plot("Both Structure factors");
	Plot variance = new Plot("Variance");
	Plot sfVsKR = new Plot("Structure Function");
	StructureFactor sf;
    //IsingLangevin ising;
    IsingField2D ising;
	Accumulator sfTheoryAcc;
	Accumulator sfTheory2Acc;
    Accumulator varianceAcc;
    Accumulator meanPhiAcc;
    
	public static void main(String[] args) {
		new Control(new LinearTheoryApp(), "Ising Field");
	}

	public LinearTheoryApp() {
		frameTogether("Plots", grid, sfVsKR, structurePeakV, variance);
		params.addm("Zoom", new ChoiceValue("Yes", "No"));
		params.addm("Interaction", new ChoiceValue("Circle", "Square"));
		params.addm("Dynamics?", new ChoiceValue("Langevin Conserve M", "Langevin No M Convervation"));
		params.add("Init Conditions", new ChoiceValue("Random Gaussian", 
				"Artificial Stripe 3", "Read From File", "Constant" ));
		params.addm("Approx", new ChoiceValue("Exact Stable",
				"Avoid Boundaries", "Exact SemiStable", "Exact", "Linear",  "Phi4"));
		params.addm("Noise", new DoubleValue(0.0, 0.0, 1.0).withSlider());
		params.addm("Horizontal Slice", new DoubleValue(0.5, 0, 0.9999).withSlider());
		params.addm("Vertical Slice", new DoubleValue(0.5, 0, 0.9999).withSlider());
		//params.addm("Conserve M?", new ChoiceValue("Yes", "No"));
		params.addm("T", 0.05);
		params.addm("H", 0.0);
		params.add("Magnetization", 0.5);
		params.addm("dt", 1.0);
		params.addm("J", -1.0);
		params.addm("R", 1000000.0);
		params.add("L/R", 4.0);
		params.add("R/dx", 16.0);
		params.add("kR bin-width", 0.1);
		params.add("Random seed", 0);
		params.add("Max Time", 350);
		params.add("Time");
		params.add("Reps");
		params.add("Mean Phi");
		params.add("Lp");
		
	}
	
	public void animate() {
		ising.readParams(params);
		//params.set("Mean Phi", ising.mean(ising.phi));
		if (params.sget("Zoom").equals("Yes")) 	grid.setAutoScale();
		else grid.setScale(-1, 1);
		structurePeakV.setAutoScale(true);
		structurePeakH.setAutoScale(true);
		sfPeakBoth.setAutoScale(true);
		sfVsKR.setAutoScale(true);
		grid.registerData(ising.Lp, ising.Lp, ising.phi);
//		if(ising.circleInt() == true){
			structurePeakV.registerLines("Peak Value", sf.getPeakC(), Color.BLACK);
			structurePeakV.registerLines("Theory", sfTheoryAcc, Color.BLUE);
			structurePeakV.registerLines("Theory2", sfTheory2Acc, Color.BLUE);
			structurePeakV.registerLines("2nd Peak", sf.get2PeakC(), Color.RED);
			variance.registerLines("Variance", varianceAcc, Color.BLACK);
			sfVsKR.registerLines("SF", sf.getAccumulatorC(), Color.BLACK);
//		}
//		}else{
//			structurePeakV.registerLines("Peak Value", sf.getPeakC(), Color.BLACK);
//			//structurePeakV.registerLines("Vertical Peak", sf.getPeakV(), Color.CYAN);
//			structurePeakV.registerLines("Theory", sfTheoryAcc, Color.BLUE);
//			structurePeakV.registerLines("Theory2", sfTheory2Acc, Color.BLUE);
//			structurePeakV.registerLines("2nd Peak", sf.get2PeakC(), Color.RED);
//			variance.registerLines("Variance", varianceAcc, Color.BLACK);
//			sfVsKR.registerLines("SF", sf.getAccumulatorC(), Color.BLACK);
//		}
	}

	public void clear() {
	}

	public void run() {
		ising = new IsingField2D(params);
		sfTheoryAcc = new Accumulator(ising.dt);
		sfTheory2Acc = new Accumulator(ising.dt);
		double binWidth = params.fget("kR bin-width");
		binWidth = IsingLangevin.KR_SP / floor(IsingLangevin.KR_SP/binWidth);
        sf = new StructureFactor(ising.Lp, ising.L, ising.R, binWidth, ising.dt);
		sf.setBounds(0.1, 14);	
		double density = findMeanPhi();
		int kR1int;
		int kR2int;
		if(params.sget("Interaction")=="Circle"){
			kR1int = getkRint(5.13562230);
			kR2int = getkRint(11.6198);
			System.out.println("circle int = " + kR1int + " " + kR2int);
		}else{
			kR1int = getkRint(4.4934092);
			kR2int = getkRint(10.9041);
			System.out.println("square int = " + kR1int + " " + kR2int);
		}
		fillTheoryAccum(density, kR1int, kR2int);
		varianceAcc = new Accumulator(params.fget("dt"));
		meanPhiAcc = new Accumulator(params.fget("dt"));
		//double kR = sf.getCircleKR();

		int reps = 0;
		while (true) {
			ising.randomizeField(density);
			for (double t = 0.0; t < params.fget("Max Time"); t = t + params.fget("dt")){
				params.set("Time", t);
				params.set("Mean Phi", ising.mean(ising.phi));
				ising.simulate();
				//accumTheoryPoint(kR, t);
				//sf.accumulateAll(t, ising.coarseGrained());
				sf.accumExact(t, ising.coarseGrained(),kR1int,kR2int);
				varianceAcc.accum(t, ising.phiVariance());
				meanPhiAcc.accum(t,ising.mean(ising.phi));
				//System.out.println(t);
				Job.animate();
			}
			reps += 1;
			params.set("Reps", reps);
		}
	}
	
//	private void accumTheoryPoint(double kR, double t) {
//		sfTheoryAcc.accum(t, linearTheory(kR, ising.lastMu, t));
//	}

	private int getkRint(double peakValue) {
		double dblePeakLength = peakValue*ising.L/(2*PI*ising.R);
		int circlePeakInt = (int)dblePeakLength;
		if(abs(2*PI*circlePeakInt*ising.R/ising.L - peakValue) >= abs(2*PI*(circlePeakInt+1)*ising.R/ising.L - peakValue))
			circlePeakInt = circlePeakInt + 1;
		return circlePeakInt;
	}

	public double findMeanPhi(){
		for (double t = 0.0; t < 50.0; t = t + params.fget("dt"))
		ising.simulate();
		return ising.mean(ising.phi);
	}
	
	public double circleLinearTheory(double kR, double mu, double time){
		//double V = ising.Lp*ising.Lp;
		double rho = ising.DENSITY;
		double D = -circlePotential(kR) - ising.T/ (1-rho*rho);
		D *= pow(1-rho*rho,2);
		//double sf = (exp(2*time*D)*(V + ising.T/D)-ising.T/D)/V;
		//double sf = (exp(2*time*D)*(V + ising.T/D)-ising.T/D)/V;
		
		double sf = exp(2*time*D);
		return sf;
	}

	public double squareLinearTheory(double kR, double mu, double time){
		//double V = ising.Lp*ising.Lp;
		double rho = ising.DENSITY;
		double D = -squarePotential(kR) - ising.T/ (1-rho*rho);
		D *= pow(1-rho*rho,2);
		//double sf = (exp(2*time*D)*(V + ising.T/D)-ising.T/D)/V;
		//double sf = (exp(2*time*D)*(V + ising.T/D)-ising.T/D)/V;
		
		double sf = exp(2*time*D);
		return sf;
	}
	
	public void fillTheoryAccum(double density, int kR1int, int kR2int){
		sfTheoryAcc = new Accumulator(ising.dt);
		sfTheory2Acc = new Accumulator(ising.dt);
		//double kR = sf.circlekRValue();
		double kR1 = ising.R*2*PI*kR1int/ising.L;
		double kR2 = ising.R*2*PI*kR2int/ising.L;
		if(params.sget("Interaction")=="Square"){
			for(double time = 0.0; time < params.fget("Max Time"); time = time + params.fget("dt")){
				sfTheoryAcc.accum(time, squareLinearTheory(kR1, 0, time));
				sfTheory2Acc.accum(time, squareLinearTheory(kR2, 0, time));
			}
		}else{
			for(double time = 0.0; time < params.fget("Max Time"); time = time + params.fget("dt")){
				sfTheoryAcc.accum(time, circleLinearTheory(kR1, 0, time));
				sfTheory2Acc.accum(time, circleLinearTheory(kR2, 0, time));
			}
		}
	}
	
	public double squarePotential(double kR){
		return sin(kR)/kR;
	}
	
	public double circlePotential(double kR){
		return 2*j1(kR)/kR;
	}

}

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
		params.addm("Interaction", new ChoiceValue("Square", "Circle"));
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
		params.add("Magnetization", 0.0);
		params.addm("dt", 0.01);
		params.addm("J", -1.0);
		params.addm("R", 1000000.0);
		params.add("L/R", 2.56);
		params.add("R/dx", 50.0);
		params.add("kR bin-width", 0.1);
		params.add("Random seed", 0);
		params.add("Max Time", 10.0);
		params.add("Time");
		params.add("Reps");
		params.add("Mean Phi");
		params.add("Lp");
		
	}
	
	public void animate() {
		ising.readParams(params);
		params.set("Time", ising.time());
		//params.set("Mean Phi", ising.mean(ising.phi));
		if (params.sget("Zoom").equals("Yes")) 	grid.setAutoScale();
		else grid.setScale(-1, 1);
		structurePeakV.setAutoScale(true);
		structurePeakH.setAutoScale(true);
		sfPeakBoth.setAutoScale(true);
		sfVsKR.setAutoScale(true);
		grid.registerData(ising.Lp, ising.Lp, ising.phi);
		if(ising.circleInt() == true){
			structurePeakV.registerLines("Peak Value", sf.getPeakC(), Color.BLACK);
			structurePeakV.registerLines("Theory", sfTheoryAcc, Color.BLUE);
			structurePeakV.registerLines("Theory2", sfTheory2Acc, Color.BLUE);
			structurePeakV.registerLines("2nd Peak", sf.get2PeakC(), Color.RED);
			variance.registerLines("Variance", varianceAcc, Color.BLACK);
			sfVsKR.registerLines("SF", sf.getAccumulatorC(), Color.BLACK);
		}else{
			structurePeakV.registerLines("Peak Value", sf.getPeakC(), Color.BLACK);
			//structurePeakV.registerLines("Vertical Peak", sf.getPeakV(), Color.CYAN);
			structurePeakV.registerLines("Theory", sfTheoryAcc, Color.BLUE);
			structurePeakV.registerLines("Theory2", sfTheory2Acc, Color.BLUE);
			structurePeakV.registerLines("2nd Peak", sf.get2PeakC(), Color.RED);
			variance.registerLines("Variance", varianceAcc, Color.BLACK);
			sfVsKR.registerLines("SF", sf.getAccumulatorC(), Color.BLACK);
		}
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
		//double density = findMeanPhi();
		double density = 0;

		//int kR1int = ising.findKRSquareInt(IsingField2D.KRsquare);//getkRint(5.13562230);
		int kR1int = ising.findKRSquareInt(5.006913291658733);
		int kR2int = ising.findKRSquareInt(IsingField2D.KRsquare2);
		double maxTime = params.fget("Max Time");
		fillTheory(density, kR1int, kR2int, maxTime);
		varianceAcc = new Accumulator(params.fget("dt"));
		meanPhiAcc = new Accumulator(params.fget("dt"));
		int reps = 0;
		double abortTimeSum = 0;
		
		while (true) {
			ising.randomizeField(density);
			ising.restartClock();
			boolean abort = false;
			collect(kR1int, kR2int);
			Job.animate();
			params.set("dt", 0.01);
			while(ising.time()< maxTime & abort == false){
				//if(ising.time()>4.5) params.set("dt", 0.001);
				params.set("Mean Phi", ising.mean(ising.phi));
				abort = ising.simulateUnstable();	
				//ising.simulate();
				if (abort){
					System.out.println("abort at " + ising.time());
					abortTimeSum += ising.time();
				}else collect(kR1int, kR2int);
				Job.animate();
			}
		reps += 1;
		params.set("Reps", reps);
		System.out.println(abort);
		double aveAbortTime = abortTimeSum/(double)reps;
		System.out.println("ave abort time = " + aveAbortTime);
		}
	}
	
	public void collect(int kR1int, int kR2int){
		sf.accumulateAll(ising.time(), ising.coarseGrained());
		sf.accumExact(ising.time(), ising.coarseGrained(),kR1int,kR2int);
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
	
	public double circleLinearTheory(double kR, double time){
		System.out.println("check this routine!!");
		//double V = ising.Lp*ising.Lp;
		double rho = ising.DENSITY;
		double D = -circlePotential(kR) - ising.T/ (1-rho*rho);
		//D *= pow(1-rho*rho,2);
		//double sf = (exp(2*time*D)*(V + ising.T/D)-ising.T/D)/V;
		//double sf = (exp(2*time*D)*(V + ising.T/D)-ising.T/D)/V;
		
		double sf = exp(2*time*D);
		return sf;
	}

	public double squareLinearTheory(double kR, double time){
		double V = ising.Lp*ising.Lp;
		double rho = ising.DENSITY;
		double D = ising.mobility*(-squarePotential(kR) - ising.T/ (1-rho*rho));
		//D *= pow(1-rho*rho,2);
		double sf = (exp(2*time*D)*(V + ising.T*ising.mobility/D)-ising.T*ising.mobility/D)/V;	
		//double sf = exp(2*time*D);
		return sf;
	}
	
	private void fillTheory(double density, int kR1int, int kR2int, double maxTime){
		sfTheoryAcc = new Accumulator(ising.dt);
		sfTheory2Acc = new Accumulator(ising.dt);
		double kR1 = ising.R*2*Math.PI*kR1int/ising.L;
		double kR2 = ising.R*2*Math.PI*kR2int/ising.L;
		ising.randomizeField(density);
		ising.restartClock();	
		double timeNotify = 0;
		if(params.sget("Interaction")=="Square"){
			sfTheoryAcc.accum(ising.time(), squareLinearTheory(kR1, ising.time()));
			sfTheory2Acc.accum(ising.time(), squareLinearTheory(kR2, ising.time()));
			while(ising.time()<maxTime){
				if(ising.time() > timeNotify){
					System.out.println("init time " + ising.time() + " of " + maxTime);
					timeNotify += 1;
				}
				ising.simulate();
				sfTheoryAcc.accum(ising.time(), squareLinearTheory(kR1, ising.time()));
				sfTheory2Acc.accum(ising.time(), squareLinearTheory(kR2, ising.time()));
			}
		}
		else if(params.sget("Interaction")=="Circle"){
			sfTheoryAcc.accum(ising.time(), circleLinearTheory(kR1, ising.time()));
			sfTheory2Acc.accum(ising.time(), circleLinearTheory(kR2, ising.time()));
			while(ising.time()<maxTime){	
				if(ising.time() > timeNotify){
					System.out.println("init time " + ising.time() + " of " + maxTime);
					timeNotify += 1;
				}
				ising.simulate();
				sfTheoryAcc.accum(ising.time(), circleLinearTheory(kR1, ising.time()));
				sfTheory2Acc.accum(ising.time(), circleLinearTheory(kR2, ising.time()));
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

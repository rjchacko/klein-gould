package rachele.demos.apps;

import static java.lang.Math.log;
import static java.lang.Math.pow;
import static scikit.util.Utilities.asList;

import java.awt.Color;

import scikit.dataset.Accumulator;
import scikit.graphics.dim2.Geom2D;
import scikit.graphics.dim2.Plot;
import scikit.jobs.Control;
import scikit.jobs.Job;
import scikit.jobs.Simulation;
import scikit.util.DoubleArray;

public class SpinodalDecomp extends Simulation{
	Plot plot = new Plot("Spinodal Decomposition");
	Accumulator freeEngAcc;
	double [] freeEnergy; 
	int noReps = 99;
	double dc = 1/(noReps+1.0);
	double concentration = .5;
	double cA = .5;
	double cB = .5;
	double T;
	
	public static void main(String[] args) {
		new Control(new SpinodalDecomp(), "Spinodal Decomposition");
	}
	
	public void load(Control c) {
		c.frame(plot);
		freeEngAcc = new Accumulator(dc);
		freeEnergy = new double[noReps];
		plot.setAutoScale(true);
		params.add("T");
	}
	

	public void animate() {
		plot.registerLines("Free Energy", freeEngAcc,  Color.black);
		double f_cA = freeEnergyCalc(T, cA);
		System.out.println(cA + " " + cB + " " + f_cA);

		plot.setDrawables(asList(
				Geom2D.line(concentration, DoubleArray.min(freeEnergy), concentration, DoubleArray.max(freeEnergy), Color.RED),
				Geom2D.line(cA, f_cA, cB, f_cA, Color.BLUE)));

	}


	public void clear() {
	}


	public void run() {
		T = 0.5;
		double dT = .01;
		int rep = 0;
		int quenchReps = 16;
		while(true){
			rep += 1;
			if(rep <= quenchReps){
				cA = .5;
				freeEnergyFill(T);
				params.set("T", T);
				Job.animate();
				T -= dT;
			}
			if (rep == quenchReps) T += dT;
			if(rep > quenchReps){
				cA = cA - dc;
				cB = 1-cA;
				Job.animate();
			}
			
		}
	}

	public void freeEnergyFill(double T){
		freeEngAcc.clear();
		for(int i = 0; i < noReps; i ++){
			double c = (i+1)*dc;
			double freeEng = freeEnergyCalc(T,c);
			//System.out.println(c + " " + freeEng);
			freeEnergy[i] = freeEng;
			freeEngAcc.accum(c,freeEng);
		}
	}
	
	public double freeEnergyCalc(double T, double c){
		double S = - c*log(c) - (1-c)*log(1-c);
		return -pow((c-.5),2)-T*S;
	}
}

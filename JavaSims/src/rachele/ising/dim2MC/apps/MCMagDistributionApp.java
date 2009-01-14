package rachele.ising.dim2MC.apps;

import static java.lang.Math.exp;
import static scikit.util.Utilities.format;

import java.awt.Color;

import rachele.ising.dim2MC.IsingLR;
import scikit.dataset.Function;
import scikit.dataset.Histogram;
import scikit.graphics.dim2.Grid;
import scikit.graphics.dim2.Plot;
import scikit.jobs.Control;
import scikit.jobs.Job;
import scikit.jobs.Simulation;
import scikit.jobs.params.ChoiceValue;
import scikit.jobs.params.DirectoryValue;

public class MCMagDistributionApp extends Simulation{
	Plot magPlot = new Plot("Magntezation");
	Grid grid = new Grid("Long Range Ising Model");
	IsingLR sim;
	int dx; 
	int accNo = 6;
	int mu; //theoretical parameter mu = 1 for Glauber and 2 for met
	double [] timeArray = new double [accNo];
	Histogram [] mag = new Histogram [accNo];	
	public static void main(String[] args) {
		new Control(new MCMagDistributionApp(), "Ferromagnet Mag Distro");
	}

	public void load(Control c) {
		c.frameTogether("stuff", grid, magPlot);
		params.add("Data Dir",new DirectoryValue("/home/erdomi/data/decision_making/java2d_mag/testruns"));
		params.addm("Dynamics", new ChoiceValue("Ising Glauber",  "Ising Metropolis"));
		params.add("Random seed", 0);
		params.add("L", 1<<7);
		params.add("R", 50);//1<<6);
		params.add("Initial magnetization", 0.0);
		params.addm("T", 0.096548444);
		params.addm("J", 1.0);
		params.addm("h", 0.0);
		params.addm("dt", 1/(double)(1<<5));
//		params.addm("dTime", 0.25);
		params.add("time");
		params.add("magnetization");
		params.add("Lp");
		params.add("Reps");
	}

	public void animate() {
		params.set("time", format(sim.time()));
		params.set("magnetization", format(sim.magnetization()));
		params.set("Lp", sim.L/dx);
		grid.registerData(sim.L/dx, sim.L/dx, sim.getField(dx));
		for (int i = 0; i < accNo; i ++){
			StringBuffer sb1 = new StringBuffer();sb1.append("Magnetization Distro"); sb1.append(i);
			float colorChunk = (float)i/(float)accNo;
			Color col = Color.getHSBColor(colorChunk, 1.0f, 1.0f);
//			PointSet normMag = new PointSet(-1,)
			magPlot.registerLines(sb1.toString(), mag[i], col);
			sb1.append("t");
			final double plotTime = timeArray [i];
			magPlot.registerLines(sb1.toString(), new Function() {
				public double eval(double m) {
					double beta = 1.0/sim.T;
					int vol = sim.L*sim.L;
					double sigmaSq = (exp(2*mu*(beta - 1.0)*plotTime)-1.0)/(vol*(beta-1.0));
					return exp(-m*m/2*sigmaSq)/(Math.sqrt(2*Math.PI*sigmaSq));
				}
			}, Color.BLACK);
		}
	}

	public void clear() {
		
	}

	public void run() {
		sim = new IsingLR(params);
		sim.randomizeField(0.0);
		for (int i =0; i < accNo; i++){mag[i] = new Histogram(0.01);}
		if (params.sget("Dynamics")=="Ising Glauber") mu = 1; else mu = 2;
		dx = 1;
		int repNo = 0;
		double dTime = params.fget("dt");
		for (int i = 0; i < accNo; i++){
			while (sim.time() < (i+1)*dTime){
				sim.step();
				Job.animate();	
			}
			timeArray[i] = sim.time();
			System.out.println("Time " + i + " = " + timeArray[i]);
		}
		while(true){
			sim.randomizeField(0.0);
			sim.restartClock();
			for (int i = 0; i < accNo; i++){
				while (sim.time() < (i+1)*dTime){
					sim.step();
					Job.animate();	
				}
				mag[i].accum(sim.magnetization());
			}
			repNo += 1;
			params.set("Reps", repNo);

		}	
	}

		
}


package rachele.ising.dim2.apps;

import static scikit.util.Utilities.format;

import java.awt.Color;

//import rachele.ising.dim2.IsingArbitraryConnect;
import rachele.ising.dim2.IsingLR;
import scikit.dataset.Accumulator;
import scikit.graphics.dim2.Grid;
import scikit.graphics.dim2.Plot;
import scikit.jobs.Control;
import scikit.jobs.Job;
import scikit.jobs.Simulation;
import scikit.jobs.params.ChoiceValue;
import scikit.jobs.params.DirectoryValue;

/**
* 
* Generic Monte Carlo Code.  Using this to build up to ferromagnetic critical exponent
* verification to build up to verification of critical exponents on various topological networks.
* 
*/

public class MCIsingApp extends Simulation{
	
	Grid grid = new Grid("Long Range Ising Model");
	Plot magPlot = new Plot("Magnetization");
	Plot mag2Plot = new Plot("Susceptibility");
	Plot chiTPlot = new Plot("Susceptibility v T");
	IsingLR sim;
//	IsingArbitraryConnect sim;
	boolean measure=false;
	int dx,n;
	double mag, aveMag, aveMagSq;
	Accumulator magAcc = new Accumulator();
	Accumulator chiAcc = new Accumulator();
	Accumulator chiTAcc = new Accumulator();
	
	public static void main(String[] args) {
		new Control(new MCIsingApp(), "Monte Carlo");
	}

	public void load(Control c) {
		c.frameTogether("MC data",grid,magPlot,mag2Plot,chiTPlot);
		params.add("Data Dir",new DirectoryValue("/home/erdomi/data/lraim/stripeToClumpInvestigation/mcResults/DO_MCdataAppBreakdown/testRuns"));
		params.addm("Dynamics", new ChoiceValue("Ising Glauber","Kawasaki Glauber", "Kawasaki Metropolis",  "Ising Metropolis"));
		params.add("Random seed", 0);
		params.add("L", 1<<7);
		params.add("R", 1);//1<<6);
		params.add("Initial magnetization", 0.0);
		params.addm("T", 1.2);
		params.addm("J", 1.0);
		params.addm("h", 0.0);
		params.addm("dt", 1.0);//1/(double)(1<<4));
		params.addm("take data",1);
		params.add("time");
		params.add("magnetization");
		params.add("Lp");
		params.add("Measuring");
		flags.add("Clear Accs");
		flags.add("Clear Aves");
		flags.add("Measure");

	}
	

	public void animate() {
		grid.setScale(-1.0, 1.0);
		grid.registerData(sim.L/dx, sim.L/dx, sim.getField(dx));
		mag = sim.magnetization();
		params.set("time", format(sim.time()));
		params.set("magnetization", format(mag));
		magAcc.accum(sim.time(),mag);
		magPlot.registerLines("Magnetization", magAcc, Color.black);
		
		aveMag = (aveMag*n+mag)/(n+1);
		aveMagSq = (aveMagSq*n+mag*mag)/(n+1);
		n+=1;
		double chi = (aveMagSq-aveMag*aveMag)/sim.T;
		chiAcc.accum(sim.time(), chi);
		
		mag2Plot.registerLines("Susceptibility", chiAcc, Color.blue);
		if(flags.contains("Clear Accs")){
			magAcc.clear();
			chiAcc.clear();
		}
		if(flags.contains("Clear Aves")){
			aveMag = 0;
			aveMagSq = 0;
			n=0;
		}
		if(flags.contains("Measure")){
			if (measure) measure=false;
			else if (measure==false) measure=true;
		}
		int measureStep = params.iget("take data");
		params.set("Measuring", measure);
		if(measure){
			if(sim.time()%measureStep==0){
				chiTAcc.accum(sim.T,chi);
				chiTPlot.registerLines("chi T", chiTAcc, Color.red);
			}
		}
		flags.clear();
		sim.setParameters(params);
	}


	public void clear() {

		
	}


	public void run() {
		chiTAcc.enableErrorBars(true);
		sim = new IsingLR(params);
//		sim = new IsingArbitraryConnect(params);
		sim.randomizeField(params.fget("Initial magnetization"));	
		dx=1;
		while(true){
			sim.step();
			Job.animate();
		}
		
	}

}

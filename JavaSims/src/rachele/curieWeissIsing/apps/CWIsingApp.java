package rachele.curieWeissIsing.apps;

import static scikit.util.Utilities.format;
import rachele.curieWeissIsing.CWIsing;
import scikit.graphics.dim2.Grid;
import scikit.jobs.Control;
import scikit.jobs.Job;
import scikit.jobs.Simulation;
import scikit.jobs.params.DirectoryValue;

/**
* 
* Generic Monte Carlo Code for Curie-Weiss Ising model
* with Gluber MC dynamics.
* 
* 
*/

public class CWIsingApp extends Simulation{

	Grid grid = new Grid("Curie-Weiss Ising Model");
	CWIsing sim;

	public static void main(String[] args) {
		new Control(new CWIsingApp(), "CW Ising Model");
	}

	public void load(Control c) {
		c.frame(grid);
		params.add("Data Dir",new DirectoryValue("/home/erdomi/data/lraim/stripeToClumpInvestigation/mcResults/DO_MCdataAppBreakdown/testRuns"));
		params.add("Random seed", 0);
		params.add("L", 1<<7);
		params.add("Initial magnetization", 0.0);
		params.addm("T", 1.2);
		params.addm("J", 1.0);
		params.addm("h", 0.0);
		params.addm("dt", 1.0);//1/(double)(1<<4));
		params.add("time");
		params.add("magnetization");

	}

	public void animate() {
		grid.setScale(-1.0, 1.0);
		grid.registerData(sim.L, sim.L, sim.getField(1));
		params.set("time", format(sim.time()));
		sim.setParameters(params);
	}


	public void clear() {
	}



	public void run() {
		sim = new CWIsing(params);
		sim.randomizeField(params.fget("Initial magnetization"));	
		while(true){
			sim.step();
			Job.animate();
		}
		
	}

}

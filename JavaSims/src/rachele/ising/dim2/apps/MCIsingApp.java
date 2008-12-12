package rachele.ising.dim2.apps;

import static scikit.util.Utilities.format;

//import rachele.ising.dim2.IsingArbitraryConnect;
import rachele.ising.dim2.IsingLR;
import scikit.graphics.dim2.Grid;
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
	IsingLR sim;
//	IsingArbitraryConnect sim;
	boolean measure=false;
	
	public static void main(String[] args) {
		new Control(new MCIsingApp(), "Monte Carlo");
	}

	public void load(Control c) {
		c.frame(grid);
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
		sim = new IsingLR(params);
//		sim = new IsingArbitraryConnect(params);
		sim.randomizeField(params.fget("Initial magnetization"));	
		while(true){
			sim.step();
			Job.animate();
		}
		
	}

}
package rachele.Networks.Apps;

import static scikit.util.Utilities.format;
import rachele.Networks.RandomGraph;
import scikit.graphics.dim2.Grid;
import scikit.jobs.Control;
import scikit.jobs.Job;
import scikit.jobs.Simulation;
//import scikit.jobs.params.DirectoryValue;

public class RandomGraphApp extends Simulation{

	Grid grid = new Grid("Random Graph");
	RandomGraph sim;

	public static void main(String[] args) {
		new Control(new RandomGraphApp(), "Random graph Ising Model");
	}

	public void load(Control c) {
		c.frame(grid);
//		params.add("Data Dir",new DirectoryValue("/home/erdomi/data/lraim/stripeToClumpInvestigation/mcResults/DO_MCdataAppBreakdown/testRuns"));
		params.add("Random seed", 1);
		params.add("L", 1<<7);
		params.add("Initial magnetization", 0.0);
		params.addm("T", 0.44444444444);
		params.addm("J", 1.0);
		params.addm("h", 0.317);
		params.addm("z", 4);
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
		sim = new RandomGraph(params);
		sim.randomizeField(params.fget("Initial magnetization"));	
		while(true){
			sim.step();
			Job.animate();
		}
		
		
	}

}

package rachele.Networks.Apps;

import static scikit.util.Utilities.format;
import rachele.Networks.SmallWorld;
import scikit.graphics.dim2.Grid;
import scikit.jobs.Control;
import scikit.jobs.Job;
import scikit.jobs.Simulation;

public class SmallWorldApp extends Simulation{
	Grid grid = new Grid("Small World");
	SmallWorld sim;

	public static void main(String[] args) {
		new Control(new SmallWorldApp(), "Small World Ising Model");
	}

	public void load(Control c) {
		c.frame(grid);
		params.add("Random seed", 1);
		params.add("L", 1<<7);
		params.add("Initial magnetization", 0.0);
		params.addm("T", 0.44444444444);
		params.addm("J", 1.0);
		params.addm("h", 0.229);
		params.addm("R", 1);
		params.addm("p",0.8);
		params.addm("dt", 1.0);
		params.add("z");
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
		sim = new SmallWorld(params);
		sim.randomizeField(params.fget("Initial magnetization"));	
		while(true){
			sim.step();
			Job.animate();
		}
		
		
	}

}

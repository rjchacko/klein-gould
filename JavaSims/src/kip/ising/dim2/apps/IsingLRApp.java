package kip.ising.dim2.apps;

import static scikit.util.Utilities.format;
import kip.ising.dim2.IsingLR;
import scikit.graphics.dim2.Grid;
import scikit.jobs.Control;
import scikit.jobs.Job;
import scikit.jobs.Simulation;
import scikit.jobs.params.ChoiceValue;


public class IsingLRApp extends Simulation {
	public static void main(String[] args) {
		new Control(new IsingLRApp(), "Ising Model");
	}
	
	Grid grid = new Grid("Coarse Grained Display");
	int dx;
	IsingLR sim;
	
	public void load(Control c) {
		c.frame(grid);
		params.addm("Dynamics", new ChoiceValue("Kawasaki Glauber", "Kawasaki Metropolis", "Ising Glauber", "Ising Metropolis"));
		params.addm("Scale colors", new ChoiceValue("False", "True"));
		params.add("Random seed", 0);
		params.add("L", 1<<8);
		params.add("R", 1<<4);
		params.add("Initial magnetization", 0.6);
		params.addm("T", 0.11);
		params.addm("J", -1.0);
		params.addm("h", 0.0);
		params.addm("dt", 0.1);
		params.add("time");
		params.add("magnetization");
	}
	
	
	public void animate() {
		params.set("time", format(sim.time()));
		params.set("magnetization", format(sim.magnetization()));
		sim.setParameters(params);
		
		if (params.sget("Scale colors").equals("False"))
			grid.setScale(-1, 1);
		else
			grid.setAutoScale();
		grid.registerData(sim.L/dx, sim.L/dx, sim.getField(dx));
	}
	
	public void clear() {
		grid.clear();
	}
	
	public void run() {
		sim = new IsingLR(params);
		sim.setField(params.fget("Initial magnetization"));
		dx = Math.max(Integer.highestOneBit(sim.R)/8, 1);
		
		double lastUpdate = 0;
		while (true) {
			while (sim.time() - lastUpdate < 2) {
				sim.step();
				Job.animate();
			}
			lastUpdate = sim.time();
			Job.animate();
		}
	}
}

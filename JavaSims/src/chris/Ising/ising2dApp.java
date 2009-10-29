package chris.Ising;

import kip.ising.dim2.Ising2D;
import scikit.graphics.dim2.Grid;
import scikit.jobs.Control;
import scikit.jobs.Job;
import scikit.jobs.Simulation;
import scikit.jobs.params.ChoiceValue;
import scikit.util.Utilities;

public class ising2dApp extends Simulation{

    Grid grid = new Grid("Ising spins");
	Ising2D model;
	double dt;
	double[] clusters;
	
	
	public static void main(String[] args) {
		new Control(new ising2dApp(), "Ising Model");
	}
	
	public void animate() {
		
		params.set("time", Utilities.format(model.time));
		model.T = params.fget("T");
		dt = params.fget("dt");
		grid.registerData(model.L1, model.L2, model.spin);
	}

	public void clear() {
		
	}

	public void load(Control c) {

		c.frame(grid);

		params.add("Seed", 77);
		params.add("Boundary", new ChoiceValue("Periodic", "Open"));
		params.add("L", 256);
		params.add("Ratio", 1.0);
		params.addm("T", 10.0);
		params.addm("dt", 10.);
		params.add("time");

	}

	public void run() {
		
		model = new Ising2D(params.iget("Seed"), params.iget("L"), (int)(params.fget("Ratio")*params.iget("L")), params.fget("T"), params.sget("Boundary").equals("Open"));
		
		while(true){
			model.step(dt);
			Job.animate();
		}
		
		//return;
	}

}

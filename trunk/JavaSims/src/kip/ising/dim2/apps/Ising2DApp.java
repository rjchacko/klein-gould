package kip.ising.dim2.apps;

import kip.ising.PercolationSite2d;
import kip.ising.dim2.Ising2D;
import scikit.graphics.dim2.Grid;
import scikit.jobs.Control;
import scikit.jobs.Job;
import scikit.jobs.Simulation;
import scikit.util.Utilities;


public class Ising2DApp extends Simulation {
    Grid grid = new Grid("Ising spins");
    Grid perc = new Grid("Perc");
	Ising2D sim;
	double dt;
	double[] clusters;
	
	public static void main(String[] args) {
		new Control(new Ising2DApp(), "Ising Model");
	}

	public void load(Control c) {
		c.frame(grid, perc);
		
		params.add("Seed", 0);
		params.add("L", 256);
		params.add("Ratio", 1.0);
		params.add("T", 0.0);
		params.addm("dt", 1.0);
		params.add("time");
		params.add("homology");
	}
	
	public void animate() {
		params.set("time", Utilities.format(sim.time));
		
		sim.T = params.fget("T");
		dt = params.fget("dt");
		grid.registerData(sim.L1, sim.L2, sim.spin);
		
		PercolationSite2d nz = new PercolationSite2d(sim.L1, sim.L2);
		nz.occupyAndBondSites(sim.spin, 1);
		for (int i = 0; i < clusters.length; i++)
			clusters[i] = nz.clusterSize(i);
		perc.registerData(sim.L1, sim.L2, clusters);
		params.set("homology", (nz.wrap_horizontal ? "horiz ":"") + (nz.wrap_vertical ? "vert ":"") + (nz.wrap_diagonal ? "diag ":""));
	}
	
	public void clear() {
		grid.clear();
		perc.clear();
	}
	
	public void run() {
		dt = params.fget("dt");
		int seed = params.iget("Seed");
		int L2 = params.iget("L");
		int L1 = (int) (L2 * params.fget("Ratio"));
		sim = new Ising2D(seed, L1, L2, params.fget("T"));	
		clusters = new double[L1*L2];
		
        while (true) {
        	sim.step(dt);            
        	Job.animate();
        }
	}
}

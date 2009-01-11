package kip.ising.dim2.apps;

import java.awt.Color;

import kip.ising.PercolationSite2d;
import kip.ising.dim2.Ising2D;
import scikit.dataset.Accumulator;
import scikit.graphics.dim2.Grid;
import scikit.graphics.dim2.Plot;
import scikit.jobs.Control;
import scikit.jobs.Job;
import scikit.jobs.Simulation;
import scikit.jobs.params.ChoiceValue;
import scikit.util.Utilities;


public class FreezingCollectorApp extends Simulation {
    Grid grid = new Grid("Ising spins");
    Plot homPlot = new Plot("Homology");
    Accumulator horz, vert, cross;
	Ising2D sim;
	int count;
	
	public static void main(String[] args) {
		new Control(new FreezingCollectorApp(), "Ising Model");
	}

	public void load(Control c) {
		c.frame(grid, homPlot);
		
		params.add("Boundary", new ChoiceValue("Open", "Periodic"));
		params.add("L", 256);
		params.add("Ratio", 1.0);
		params.add("Count max", 10000);
		params.add("Time max", 100);
		params.add("count");
		params.add("time");
	}
	
	public void animate() {
		params.set("count", count);
		params.set("time", Utilities.format(sim.time));
		
		grid.registerData(sim.L1, sim.L2, sim.spin);
		homPlot.registerPoints("Horiz", horz, Color.BLUE);
		homPlot.registerPoints("Vert", vert, Color.RED);
		homPlot.registerPoints("Cross", cross, Color.BLACK);
	}
	
	public void clear() {
		grid.clear();
		homPlot.clear();
	}
	
	public void run() {
		boolean openBoundary = params.sget("Boundary").equals("Open");
		int L2 = params.iget("L");
		int L1 = (int) (L2 * params.fget("Ratio"));
		
		int cntMax = params.iget("Count max");
		int timeMax = params.iget("Time max");
		int sampleInterval = 10; 
		
		horz = new Accumulator();
		vert = new Accumulator();
		cross = new Accumulator();
		horz.enableErrorBars(true);
		vert.enableErrorBars(true);
		cross.enableErrorBars(true);
		
		for (count = 0; count < cntMax; count++) {
			sim = new Ising2D(count, L1, L2, 0);
			sim.openBoundary = openBoundary;
			
			while (sim.time < timeMax) {
				sim.step(1);            
				Job.animate();
				
				if (sim.time % sampleInterval == 0) {
		    		PercolationSite2d nz = new PercolationSite2d(sim.L1, sim.L2, openBoundary);
		    		nz.occupyAndBondSites(sim.spin, 1);
		    		nz.findHomologies();
		    		horz.accum(sim.time, nz.horizontalHomology() ? 1 : 0);
		    		vert.accum(sim.time, nz.verticalHomology() ? 1 : 0);
		    		cross.accum(sim.time, (nz.crossHomology() || nz.pointHomology()) ? 1 : 0);
				}
 			}
		}
	}
}

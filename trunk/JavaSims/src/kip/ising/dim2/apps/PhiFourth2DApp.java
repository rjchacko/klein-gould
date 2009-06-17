package kip.ising.dim2.apps;

import static scikit.util.Utilities.format;
import kip.ising.dim2.PhiFourth2D;
import scikit.graphics.dim2.Grid;
import scikit.jobs.Control;
import scikit.jobs.Job;
import scikit.jobs.Simulation;
import scikit.jobs.params.ChoiceValue;
import scikit.jobs.params.DoubleValue;

public class PhiFourth2DApp extends Simulation {
	Grid grid = new Grid("Grid");
	Grid gmode = new Grid("Growth Mode");
//	Plot plot = new Plot("Slice");
	PhiFourth2D clump;
	
	public static void main(String[] args) {
		new Control(new PhiFourth2DApp(), "Clump Model Saddle Profile");
	}

	public void load(Control c) {
		c.frame(grid);
		params.addm("Saddle", new ChoiceValue("No", "Yes"));
		params.addm("Noise", new ChoiceValue("Yes", "No"));
		params.addm("T", new DoubleValue(0.05, -0.2, 0.1).withSlider());
		params.addm("h", new DoubleValue(0., -0.2, 0.2).withSlider());
		params.addm("dt", 1.0);
		params.add("R", 1000.0);
		params.add("L", 8000.0);
		params.add("dx", 100.0);
		params.add("Random seed", 0);
		params.add("Time");
		params.add("F density");
		params.add("dF/dphi");
		params.add("Rx");
		params.add("Ry");
		params.add("Eigenvalue");
		params.add("del Eigenmode");
		params.add("Valid profile");
		flags.add("Res up");
		flags.add("Res down");
	}
	
	public void animate() {
		if (flags.contains("Res up"))
			clump.doubleResolution();
		if (flags.contains("Res down"))
			clump.halveResolution();
		flags.clear();
		clump.readParams(params);
		
		grid.setAutoScale();
		grid.setDrawRange(true);
		int Lp = clump.numColumns();
		grid.registerData(Lp, Lp, clump.phi());
		
//		gmode.registerData(Lp, Lp, clump.growthEigenmode);
//		params.set("Eigenvalue", format(clump.growthEigenvalue));
//		params.set("del Eigenmode", format(clump.rms_growthEigenmode));
		
//		double[] section = new double[Lp];
//		System.arraycopy(clump.growthEigenmode, Lp*(Lp/2), section, 0, Lp);
//		plot.registerLines("", new PointSet(0, 1, section), Color.BLUE);
		
		params.set("dx", clump.dx);
		params.set("Time", format(clump.time()));
		params.set("F density", format(clump.freeEnergyDensity));
		params.set("dF/dphi", format(clump.rms_dF_dphi));
		params.set("Rx", clump.Rx);
		params.set("Ry", clump.Ry);
		params.set("Valid profile", !clump.rescaleClipped);
	}
	
	public void clear() {
		grid.clear();
//		plot.clear();
	}
	
	public void run() {
		clump = new PhiFourth2D(params);
		clump.initializeFieldWithRandomSeed();
		Job.animate();
		
		while (true) {			
			if (params.sget("Saddle").equals("Yes")) {
				clump.saddleStep();
			}
			else {
				clump.simulate();
				clump.relaxInteraction();
			}
			
			Job.animate();
		}
	}

}
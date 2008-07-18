package rachele.ising.dim2.apps;

import static scikit.util.Utilities.format;

import java.awt.Color;

import rachele.ising.dim2.IsingLR;
import rachele.ising.dim2.StructureFactor;
import scikit.graphics.dim2.Grid;
import scikit.graphics.dim2.Plot;
import scikit.jobs.Control;
import scikit.jobs.Job;
import scikit.jobs.Simulation;
import scikit.jobs.params.ChoiceValue;

public class IsingLRStructApp extends Simulation {
	public static void main(String[] args) {
		new Control(new IsingLRStructApp(), "Ising Model");
	}
	
	Grid fieldDisplay = new Grid("Coarse Grained Display");
	//FieldDisplay fieldDisplay = new FieldDisplay("Coarse Grained Display", true);
	Plot structureDisplayH = new Plot("Structure Factor - Vertical Component");
	Plot structureDisplayV = new Plot("Structure Factor - Horizontal Component");
	Plot circleStructureDisplay = new Plot("Structure Factor - Circle Average");
    //Plot hSlice = new Plot("Horizontal Slice", true);
    //Plot vSlice = new Plot("Vertical Slice", true);

	int dx;
	StructureFactor structure;
	IsingLR sim;
	
	public void load(Control c) {
		c.frameTogether("Plots", fieldDisplay, structureDisplayH, structureDisplayV, circleStructureDisplay);
		params.addm("Dynamics", new ChoiceValue("Kawasaki Glauber", "Kawasaki Metropolis", "Ising Glauber", "Ising Metropolis"));
		params.addm("Scale colors", new ChoiceValue("False", "True"));
		//params.addm("Horizontal Slice", new DoubleValue(0.5, 0, 0.9999).withSlider());
		//params.addm("Vertical Slice", new DoubleValue(0.5, 0, 0.9999).withSlider());
		params.add("Random seed", 0);
		params.add("L", 1<<8);
		params.add("R", 1<<4);
		params.add("Initial magnetization", 0.6);
		params.addm("T", 0.11);
		params.addm("J", -1.0);
		params.addm("h", 0.0);
		params.addm("dt", 1.0);
		params.addm("init time", 0);
		params.add("time");
		params.add("magnetization");
		params.add("Lp");
		
		flags.add("Clear S.F.");
	}
	
	
	public void animate() {
		params.set("time", format(sim.time()));
		params.set("magnetization", format(sim.magnetization()));
		params.set("Lp", sim.L/dx);
		sim.setParameters(params);
	
		
		//fieldDisplay.setData(sim.L/dx, sim.L/dx, sim.getField(dx));
		fieldDisplay.registerData(sim.L/dx, sim.L/dx, sim.getField(dx));
		//if (params.sget("Scale colors").equals("False"))
		//	fieldDisplay.setScale(-1, 1);
		//else
		//	fieldDisplay.setAutoScale();
		structureDisplayV.registerLines("vertical", structure.getAccumulatorV(), Color.BLACK);
		structureDisplayH.registerLines("horizontal", structure.getAccumulatorH(), Color.BLACK);
		structureDisplayV.registerLines("vertical ave", structure.getAccumulatorVA(), Color.BLUE);	
		structureDisplayH.registerLines("horizontal ave", structure.getAccumulatorHA(), Color.BLUE);
		circleStructureDisplay.registerLines("circle", structure.getAccumulatorC(), Color.RED);	
		circleStructureDisplay.registerLines("circle ave", structure.getAccumulatorCA(), Color.YELLOW);	

	
		if (flags.contains("Clear S.F.")) {
			structure.getAccumulatorCA().clear();
			structure.getAccumulatorHA().clear();
			structure.getAccumulatorVA().clear();
			System.out.println("clicked");
		}
		flags.clear();
	}
	
	
	public void clear() {
	}
	
	
	public void run() {
		sim = new IsingLR(params);
		sim.setField(params.fget("Initial magnetization"));
		dx = Math.max(Integer.highestOneBit(sim.R)/8, 1);
		structure = new StructureFactor(sim.L/dx, sim.L, sim.R, 0.1, sim.dTime());
		
		System.out.println("equilibrating");
		//while (sim.time() < params.fget("init time")) {
		//	sim.step();
		//	Job.animate();
		//}
		
		System.out.println("running");
		//double lastUpdate = 0;
		while (true) {
			//while (sim.time() - lastUpdate < 2) {
				sim.step();
				//Job.animate();
			//}
			//lastUpdate = sim.time();
			structure.getAccumulatorH().clear();
			structure.getAccumulatorV().clear();			
			structure.getAccumulatorC().clear();
			structure.accumulateAll(sim.time(), sim.getField(dx));
			//avStructH.accum(structure.getAccumulatorH());
			Job.animate();
		}
	}
}

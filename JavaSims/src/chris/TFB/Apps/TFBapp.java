package chris.TFB.Apps;

import java.awt.Color;
import java.text.DecimalFormat;

import scikit.graphics.ColorPalette;
import scikit.graphics.dim2.Grid;
import scikit.graphics.dim2.Plot;
import scikit.jobs.Control;
import scikit.jobs.Job;
import scikit.jobs.Simulation;
import scikit.jobs.params.ChoiceValue;
import scikit.jobs.params.DirectoryValue;
import chris.TFB.TFBmodel;

public class TFBapp extends Simulation{
	
	private TFBmodel model;
	private boolean draw;
	private int L;
	
	private Grid grid = new Grid("Fibre State");
	private ColorPalette palette1;
	
	private DecimalFormat fmt = new DecimalFormat("0.000");
	
	private Plot gplot = new Plot("Free Energy");
	
	/////////////////////////////////////////////////////////
	
	public static void main(String[] args) {
		new Control(new TFBapp(), "TFB Model");
	}
	
	
	public void load(Control c) {
		params.add("Data Directory",new DirectoryValue("/Users/cserino/Desktop/"));
		params.add("Random Seed",0);
		params.add("Lattice Size",1<<8);
		params.add("Stress / site", 0.1);		
		params.add("Failure Stress",5.0);
		params.add("\u03C3_f width",(double) 0);
		params.add("Kappa", (double) 0.5);
		params.add("D", (double) 0.5);
		params.add("Temperature", (double) 0.2);
		params.add("Animation", new ChoiceValue("Off","On"));
		params.addm("Record", new ChoiceValue("Off","On"));
		params.add("MC Time Step");
		params.add("Phi");
		
		c.frameTogether("TFB",grid,gplot);
	}
	
	public void run() {
	
		// initialize model
		model = new TFBmodel(params);
		
		// setup grid
		setupColors();
		
		// Run the simulation
		Job.animate();
		while(model.nextLattice()){
			model.takeData();
			Job.animate();
		}
		
	}
	
	public void setupColors(){
		
		model.writeHeader();
		L = model.getL();
		draw = (params.sget("Animation").equals("On"));
		palette1  = new ColorPalette();
		palette1.setColor(0,Color.BLACK);
		palette1.setColor(1,Color.WHITE);
		grid.setColors(palette1);
		
		return;
	}
	
	public void animate() {
		
		params.set("MC Time Step", model.getTime());
		params.set("Phi", fmt.format(model.getPhi()));

		if(!draw) return;
		
		grid.registerData(L,L,model.getState());
		gplot.registerLines("Free Energy",model.getGplot(), Color.BLACK);

		if(params.sget("Record").equals("On")) model.takePicture(grid);

	}

	
	public void clear() {
		
		grid.clear();
		return;
	}
	
	
}
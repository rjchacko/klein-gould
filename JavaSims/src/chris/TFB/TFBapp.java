package chris.TFB;

import java.awt.Color;
import java.text.DecimalFormat;

import scikit.graphics.ColorPalette;
import scikit.graphics.dim2.Grid;
import scikit.jobs.Control;
import scikit.jobs.Job;
import scikit.jobs.Simulation;
import scikit.jobs.params.ChoiceValue;
import scikit.jobs.params.DirectoryValue;
import scikit.jobs.params.DoubleValue;

public class TFBapp extends Simulation{
	
	private TFBmodel model;
	private boolean draw;
	private int L;
	
	private Grid grid = new Grid("Fibre State");
	private ColorPalette palette1;
	
	private DecimalFormat fmt = new DecimalFormat("0.000");
	
	/////////////////////////////////////////////////////////
	
	public static void main(String[] args) {
		new Control(new TFBapp(), "TFB Model");
	}
	
	
	public void load(Control c) {
		params.add("Data Directory",new DirectoryValue("/Users/cserino/Desktop/"));
		params.add("Random Seed",0);
		params.add("Lattice Size",1<<8);
		params.add("Stress / site",1.0);		
		params.add("Failure Stress",5.0);
		params.add("\u03C3_f width",(double) 0);
		params.add("Kappa", (double) 1);
		params.addm("Temperature", new DoubleValue(1, 0, 5000).withSlider());
		params.add("Animation", new ChoiceValue("Off","On"));
		params.addm("Record", new ChoiceValue("Off","On"));
		params.add("MC Time Step");
		params.add("Phi");
		
		c.frame(grid);
	}
	
	public void run() {
	
		// initialize model
		model = new TFBmodel(params);
		
		// setup grid
		setupColors();
		
		// Run the simulation
		while(model.nextLattice()){
			model.takeData();
			Job.animate();
		}
		
	}
	
	public void setupColors(){
		
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
		if(params.sget("Record").equals("On")) model.takePicture(grid);

	}

	
	public void clear() {
		
		grid.clear();
		return;
	}
	
	
}

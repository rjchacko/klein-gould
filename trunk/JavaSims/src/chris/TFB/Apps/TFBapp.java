package chris.TFB.Apps;

import java.awt.Color;
import java.io.File;
import java.text.DecimalFormat;

import scikit.dataset.PointSet;
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
	private Plot gplot = new Plot("Free Energy");
	private PointSet gEnergy;
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
		params.add("Stress / site", 1E-5);		
		params.add("Failure Stress",5.0);
		params.add("\u03C3_f width",(double) 0);
		params.add("Kappa", (double) 0.5);
		params.add("D", (double) 0.5);
		params.add("Temperature", (double) 100);
		params.add("Animation", new ChoiceValue("On","Off"));
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
		
		// draw the free energy
		drawFE();
		
		// write parameters to file
		TFBmodel.printParams(params.sget("Data Directory")+File.separator+"Params.txt", params);
		
		// Run the simulation
		Job.animate();
		while(model.nextLattice()){
			model.takeData();
			Job.animate();
		}
		
	}
	
	public void drawFE(){
		
		double[] crd = new double[model.getN()+1];
		for(int jj = 0 ; jj < model.getN() + 1 ; jj++){
			crd[jj] = ((double)(jj))/((double) model.getN());
		}
		
		//gEnergy = new PointSet(crd,model.gLandscape());
		gEnergy = new PointSet(crd,model.hLandscape());

		
		return;
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
		gplot.registerLines("Free Energy",gEnergy, Color.BLACK);

		if(params.sget("Record").equals("On")) model.takePicture(grid);

	}

	
	public void clear() {
		
		grid.clear();
		return;
	}
	
	
}

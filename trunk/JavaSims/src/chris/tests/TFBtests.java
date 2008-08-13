package chris.tests;

import java.awt.Color;

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
import chris.util.PrintUtil;

public class TFBtests extends Simulation{
	
	private TFBmodel model;
	private boolean draw;
	private int L;
	
	private Grid grid = new Grid("Fibre State");
	private ColorPalette palette1;
	
//	private DecimalFormat fmt = new DecimalFormat("0.000");
	
	private Plot gplot = new Plot("Free Energy");
	
	PointSet gEnergy;

	
	/////////////////////////////////////////////////////////
	
	public static void main(String[] args) {
		new Control(new TFBtests(), "TFB Model");
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
		
		// testing free energy
		double[] crd = new double[model.getN()+1];
		double[] tmp = new double[model.getN()+1];	
		
		for(int jj = 0 ; jj < model.getN()+1 ; jj++){
			crd[jj] = jj/model.getN();
			tmp[jj] = model.gEnergy(crd[jj]);
			//tmp[jj] = 1E-5*crd[jj];
		}
		
		PrintUtil.printArrayToFile("/Users/cserino/Desktop/debugX.txt",crd,model.getN()+1,1);
		PrintUtil.printArrayToFile("/Users/cserino/Desktop/debugY.txt",tmp,model.getN()+1,1);
		
		gEnergy = new PointSet(crd,tmp);
		
		Job.animate();
		
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
		
//		params.set("MC Time Step", model.getTime());
//		params.set("Phi", fmt.format(model.getPhi()));

		if(!draw) return;
		
		grid.registerData(L,L,model.getState());
		gplot.registerLines("Free Energy",gEnergy, Color.BLACK);

		//if(params.sget("Record").equals("On")) model.takePicture(grid);

	}

	
	public void clear() {
		
		grid.clear();
		return;
	}
	
	
}

package chris.ofcdamage;

import java.awt.Color;
import java.io.File;
import java.text.DecimalFormat;

import scikit.graphics.ColorGradient;
import scikit.graphics.ColorPalette;
import scikit.graphics.dim2.Grid;
import scikit.jobs.Control;
import scikit.jobs.Job;
import scikit.jobs.Simulation;
import scikit.jobs.params.ChoiceValue;
import scikit.jobs.params.DirectoryValue;
import scikit.jobs.params.DoubleValue;
import chris.ofc.Damage2D;

public class varApp extends Simulation{

	varLives model;
	
	Grid gridS = new Grid ("Stress");
	Grid gridL = new Grid ("Lives");
	
	ColorPalette palette1;
	ColorGradient cGradient;
	
	DecimalFormat fmt = new DecimalFormat("0000.00");
	
	private double ScMax;
	private Boolean pretime;
	private int pt, ptmax;
	private boolean draw;
	
	public static void main(String[] args) {
		new Control(new varApp(), "OFC Parameters");
	}
	
public void load(Control c) {
		
		params.add("Data Directory",new DirectoryValue("/Users/cserino/Desktop/"));
		params.add("Random Seed",0);
		params.add("Interaction Shape", new ChoiceValue("Circle","Square","Diamond","All Sites"));
		params.add("Interaction Radius (R)",(int)(30));
		params.add("Lattice Size",1<<8);
		params.add("Boundary Condtions", new ChoiceValue("Periodic","Bordered"));
		params.add("Min Lives", 2);	
		params.add("Max Lives", 10);
		params.add("Failure Stress (\u03C3_f)",2.0);
		params.add("\u03C3_f width",(double)(0));
		params.add("Residual Stress (\u03C3_r)",1.25);
		params.add("\u03C3_r width",0.5);
		params.add("Dissipation (\u03B1)",new DoubleValue(0.01,0,1));
		params.add("\u03B1 Width", (double)(0));
		params.add("Animation", new ChoiceValue("Off","On"));
		params.addm("Record", new ChoiceValue("Off","On"));
		params.add("Number of Plate Updates");
		params.add("Last Avalanche Size");
		params.add("N_dead");
		
		c.frameTogether("OFC Model with Damage", gridS, gridL);
		
	}
	
	public void run() {
		
		// Setup model
		params.set("N_dead","Initializing");
		model = new varLives(params);
		
		// Setup display
		setupDisplays();
		
		// Setup output files
		Damage2D.PrintParams(model.getOutdir()+File.separator+"Params.txt",params);	
		model.writeDataHeaders();
		
		// Run the simulation
		pretime = true;
		pt    = 0;
		ptmax = 100000;
		
		// Equilibrate
		model.setEquil(true);
		while(pt < ptmax){
			model.equilibrate();
			pt++;
			Job.animate();
		}
		
		// Simulate for Data
		pretime = false;
		model.setEquil(false);
		while(model.avalanche()){
			model.takeData();
			Job.animate();
			if (params.sget("Record").equals("On")){
				// FIX ME!!!!!!!!!!!!
				model.takePicture(gridS, true);
				model.takePicture(gridL, true);
			}
		}
		
		params.set("N_Dead","Finished");
		Job.animate();
		
		return;
	}
	
	private void setupDisplays(){
		
		if(draw = (params.sget("Animation").equals("On"))){
		
			// Setup color scheme
			palette1  = new ColorPalette();
			cGradient = new ColorGradient();

			if(params.sget("Mode").equals("Damage")){
				int MostLives = model.maxLives();
				Color[] Carray = new Color[]{Color.YELLOW,Color.RED,Color.GREEN,Color.BLUE,Color.GRAY};		
				palette1.setColor(0,Color.BLACK);
				for (int jj = 1 ; jj <= MostLives ; jj++){
					palette1.setColor(jj,Carray[jj%5]);
				}

				gridL.setColors(palette1);
			}
		}
		
		return;
	}

	public void animate() {

		if(pretime){
			params.set("Number of Plate Updates", pt - ptmax);
			params.set("Last Avalanche Size",model.getAvlnchSize());
			params.set("N_dead", "Equilibrating");
		}
		else{
			params.set("Number of Plate Updates",model.getTime(1));
			params.set("Last Avalanche Size", model.getAvlnchSize());
			params.set("N_dead", model.getNdead());
		}
		
		if(draw){
			
			
			int L = model.getL();
			int N = L*L;
			
			int[] foo = model.getLives();

			double[] copyStress = model.getStress();
			
			for (int jj=0 ; jj < N ; jj++){
				cGradient.getColor(copyStress[jj],-2,ScMax);
			}

			gridS.setColors(cGradient);
			gridS.registerData(L,L,copyStress);
			gridL.registerData(L, L, foo);
			
		}
	
	}

	public void clear() {
		gridS.clear();
		gridL.clear();
	}


}

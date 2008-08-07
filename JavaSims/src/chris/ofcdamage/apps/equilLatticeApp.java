package chris.ofcdamage.apps;

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
import chris.foo.ofc.Damage2D;
import chris.ofcdamage.damage;

public class equilLatticeApp extends Simulation{

	damage model;
	
	Grid gridS = new Grid ("Stress");
	Grid gridL = new Grid ("Lives");
	
	ColorPalette palette1;
	ColorGradient cGradient;
	
	DecimalFormat fmt  = new DecimalFormat("0000.00");
	DecimalFormat ifmt = new DecimalFormat("000");
	
	private double ScMax;
	private Boolean pretime;
	private int pt, ptmax;
	private boolean draw;
	
	public static void main(String[] args) {
		new Control(new equilLatticeApp(), "OFC Parameters");
	}
	
public void load(Control c) {
		
		params.add("Data Directory",new DirectoryValue("/Users/cserino/Desktop/"));
		params.add("Random Seed",0);
		params.add("Interaction Shape", new ChoiceValue("Circle","Square","Diamond","All Sites"));
		params.add("Interaction Radius (R)",(int)(30));
		params.add("Lattice Size",1<<8);
		params.add("Boundary Condtions", new ChoiceValue("Periodic","Bordered"));
		params.add("Intitial Stess", new ChoiceValue("Random","Constant"));
		params.add("Min Lives", 2);	
		params.add("Max Lives", 10);
		params.add("Failure Stress (\u03C3_f)",2.0);
		params.add("\u03C3_f width",(double)(0));
		params.add("Residual Stress (\u03C3_r)",1.25);
		params.add("\u03C3_r width",0.5);
		params.add("Dissipation (\u03B1)",new DoubleValue(0.01,0,1));
		params.add("\u03B1 width", (double)(0));
		params.add("Equil Time", 100000);
		params.add("Number to Store", 100000);
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
		model = new damage(params);
		
		// Setup display
		setupDisplays();
		
		// Setup output files
		configOF();
		
		// Run the simulation
		
		// Equilibrate
		equil();
		// Save equilibrated lattices to file 
		createEquil();

		// Fin!
		params.set("N_dead","Finished");
		Job.animate();
		
		return;
	}
	
	private void configOF(){
		
		Damage2D.PrintParams(model.getOutdir()+File.separator+"Params.txt",params);	
		model.writeDataHeaders();
		
		return;
	}
	
	private void equil(){
	
		// Equilibrate
		pretime = true;
		pt      = 0;
		ptmax   = params.iget("Equil Time");
		
		model.setEquil(true);
		model.runClocks(false);
		while(pt < ptmax){
			model.equilibrate();
			pt++;
			Job.animate();
		}
		
		return;
	}
	
	private void createEquil(){
		
		// Simulate w/o Damage for Data
		int elcounter = 0;
		params.set("N_dead","Storing Lattice(s)");
		model.setEquil(true);
		model.runClocks(true);		
		int pt0 = 0;
		ptmax = 500*params.iget("Number to Store"); 
		pt    = 500*params.iget("Number to Store");
		while(pt0 < ptmax){
			model.avalanche();
			pt++;
			if((pt0 % 500) == 0){
				// print all data to file
				model.printStress(model.getOutdir()+File.separator+"Stress_"+ ifmt.format(elcounter) +".txt");
				model.printSr(model.getOutdir()+File.separator+"Sr"+ ifmt.format(elcounter) +".txt");
				elcounter++;
			}
			pt0++;
			Job.animate();
		}
		
		return;
	}
	
	private void setupDisplays(){
		
		if(draw = (params.sget("Animation").equals("On"))){
		
			// Setup color scheme
			palette1  = new ColorPalette();
			cGradient = new ColorGradient();

			int MostLives = model.maxLives();
			Color[] Carray = new Color[]{Color.YELLOW,Color.RED,Color.GREEN,Color.BLUE,Color.GRAY};		
			palette1.setColor(0,Color.BLACK);
			for (int jj = 1 ; jj <= MostLives ; jj++){
				palette1.setColor(jj,Carray[jj%5]);
			}

			gridL.setColors(palette1);
			
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

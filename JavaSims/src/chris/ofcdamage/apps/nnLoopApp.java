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
import chris.ofcdamage.damageCluster;
import chris.util.PrintUtil;

public class nnLoopApp extends Simulation{

	damageCluster model;
	
	Grid gridS = new Grid ("Stress");
	Grid gridL = new Grid ("Lives");
	Grid gridC = new Grid ("Clusters");
	
	ColorPalette palette1, palette2;
	ColorGradient cGradient;
		
	DecimalFormat fmt = new DecimalFormat("0000.00");
	
	private double ScMax;
	private Boolean pretime;
	private int pt, ptmax, cycle;
	private boolean draw;
	
	public static void main(String[] args) {
		new Control(new nnLoopApp(), "OFC Parameters");
	}
	
public void load(Control c) {
		
		params.add("Data Directory",new DirectoryValue("/Users/cserino/Desktop/"));
		params.add("Random Seed",0);
		params.add("Interaction Shape", new ChoiceValue("Diamond","Circle","Square","All Sites"));
		params.add("Interaction Radius (R)",(int)(1));
		params.add("Lattice Size",1<<8);
		params.add("Boundary Condtions", new ChoiceValue("Periodic","Bordered"));
		params.add("Intitial Stess", new ChoiceValue("Random","Constant"));
		params.add("Min Lives", 1);	
		params.add("Max Lives", 1);
		params.add("Failure Stress (\u03C3_f)",2.0);
		params.add("\u03C3_f width",(double)(0));
		params.add("Residual Stress (\u03C3_r)",1.25);
		params.add("\u03C3_r width",0.5);
		params.add("Dissipation (\u03B1)",new DoubleValue(0.35,0,1));
		params.add("\u03B1 width", (double)(0));
		params.add("Equil Time", 1000);
		params.add("Trend Time", 0);
		params.add("Number of Cycles", (int) 1);
		params.add("Animation", new ChoiceValue("Off","On"));
		params.addm("Record", new ChoiceValue("Off","On"));
		params.add("Number of Plate Updates");
		params.add("Last Avalanche Size");
		params.add("N_dead");
		params.add("Cycles Left");
		
		c.frameTogether("OFC Model with Damage", gridS, gridL, gridC);
		
	}
	
	public void run() {
		
		cycle = params.iget("Number of Cycles");
		
		while(cycle > 0){
			cycle--;
			params.set("Cycles Left",cycle);
			//////////
			
			// Setup model
			params.set("N_dead","Initializing");
			model    = new damageCluster(params);
			
			// Setup display
			setupDisplays();
			
			// Setup output files
			configOF();
			
			// Run the simulation
			
			// Equilibrate
			equil();
			// Simulate w/o Damage for Data
			ideal();
			// Simulate w/ Damage for Data
			damage();
			
			// Fin!
			params.set("N_dead","Finished");
			Job.animate();
			
			// Save Cluster Grid for Analysis
			recordClusters();
			
			///////////
			
			params.set("Random Seed", params.iget("Random Seed") + 1);
		}
		
		return;
	}
	
	private void configOF(){
		
		Damage2D.PrintParams(model.getOutdir()+File.separator+"Params_"+(cycle+1)+".txt",params);	
		model.setDataFile(cycle+1);
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
	
	private void ideal(){
		
		// Simulate w/o Damage for Data
		params.set("N_dead","Est. Equilib");
		model.setEquil(true);
		model.runClocks(true);		
		int pt0 = 0;
		ptmax = params.iget("Trend Time"); 
		pt    = params.iget("Trend Time");
		while(pt0 < ptmax){
			model.avalanche();
			model.takeData();
			pt++;
			pt0++;
			Job.animate();
		}
		
		return;
	}
	
	private void damage(){
		
		Job.animate();
		pretime = false;
		model.setEquil(false);
		model.runClocks(true);
		while(model.avalanche()){
			model.takeData();
			Job.animate();
			if (params.sget("Record").equals("On")){
				// FIX ME!!!!!!!!!!!!  WHY DOES THIS NEED FIXING??????
				model.takePicture(gridS, true);
				model.takePicture(gridL, true);
				model.takePicture(gridC, true);
			}
		}
		
		return;
	}
	
	private void setupDisplays(){
		
		if(draw = (params.sget("Animation").equals("On"))){
		
			// Setup color scheme
			palette1  = new ColorPalette();
			palette2  = new ColorPalette();
			cGradient = new ColorGradient();

			int MostLives = model.maxLives();
			Color[] Carray = new Color[]{Color.YELLOW,Color.RED,Color.GREEN,Color.BLUE,Color.GRAY};		
			palette1.setColor(0,Color.BLACK);
			for (int jj = 1 ; jj <= MostLives ; jj++){
				palette1.setColor(jj,Carray[jj%5]);
			}

			gridL.setColors(palette1);
		
			Carray = new Color[]{Color.RED, Color.ORANGE, Color.YELLOW, Color.GREEN, Color.BLUE,
								 Color.GRAY, Color.PINK, Color.MAGENTA, Color.CYAN, Color.PINK};    
		    palette2.setColor(0,Color.WHITE);
		    for (int jj = 1 ; jj <= model.getL()*model.getL() ; jj++){
		    	palette2.setColor(jj,Carray[jj%10]);
		    }
		    
		    gridC.setColors(palette2);
			
			
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
						
			double[] copyStress = model.getStress();
			for (int jj=0 ; jj < N ; jj++){
				cGradient.getColor(copyStress[jj],-2,ScMax);
			}

			gridS.setColors(cGradient);
			gridS.registerData(L,L,copyStress);
		
			if(!pretime){
				int[] foo  = model.getLives();
				int[] foo2 = model.getClusters();
				gridL.registerData(L, L, foo);
				gridC.registerData(L, L, foo2);
			}
			
		}
	
	}

	public void clear() {
		gridS.clear();
		gridL.clear();
		gridC.clear();
	}
	
	private void recordClusters(){
		
		PrintUtil.printlnToFile(model.getOutdir()+File.separator+"Clusters_"+(cycle+1)+".txt", 
								"Percolating Cluster Number", model.pcnPC());
		PrintUtil.printVectorToFile(model.getOutdir()+File.separator+"Clusters_"+(cycle+1)+".txt",
								   model.getClusters());
		model.printCdistance(cycle);
		if(draw) model.takePicture(gridC, true, model.getOutdir()+File.separator+"PercolatingClusters_"+(cycle+1)+".png");
		return;
	}


}

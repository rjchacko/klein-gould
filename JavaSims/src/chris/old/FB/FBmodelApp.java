package chris.old.FB;

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
import chris.old.ofc.Damage2D;

public class FBmodelApp extends Simulation{

	Damage2D model;
	
	Grid gridS = new Grid ("Stress");
	Grid gridL = new Grid ("Lives");
	
	ColorPalette palette1;
	ColorGradient cGradient;
	
	DecimalFormat fmt = new DecimalFormat("0000.00");
	
	double ScMax;
	Boolean pretime;
	
	public static void main(String[] args) {
		new Control(new FBmodelApp(), "OFC Parameters");
	}
	
public void load(Control c) {
		
		params.add("Data Directory",new DirectoryValue("/Users/cserino/Desktop/"));
		params.add("Mode");
		params.add("Random Seed",0);
		params.add("Interaction Shape", new ChoiceValue("All Sites","Circle","Square","Diamond"));
		params.add("Interaction Radius (R)",(int)(30));
		params.add("Lattice Size",1<<8);
		params.add("Boundary Condtions", new ChoiceValue("Periodic","Bordered"));
		params.add("Equilibrate Time");
		params.add("Number of Lives");
		params.add("Life Distribution", new ChoiceValue("Constant","Gaussian", "Step Down"));
		params.add("LD Width",0.0);
		params.add("Failure Stress (\u03C3_f)",1.0);
		params.add("\u03C3_f width",(double)(0));
		params.add("Residual Stress (\u03C3_r)",0);
		params.add("\u03C3_r width",0);
		params.add("Dissipation (\u03B1)",new DoubleValue(.1,0,1));
		params.add("\u03B1 Width", (double)(0));
		params.add("Animation", new ChoiceValue("Off","On"));
		params.addm("Record", new ChoiceValue("Off","On"));
		params.add("Number of Plate Updates");
		params.add("Last Avalanche Size");
		params.add("<Lives Left>");
		
		params.set("Mode","Damage");
		params.set("Equilibrate Time",0);
		params.set("Number of Lives",1);
		
		c.frameTogether("OFC Model with Damage", gridS, gridL);
		
	}
	
	public void animate() {


		params.set("Number of Plate Updates",model.getTime(1));
		params.set("Last Avalanche Size",model.getAS());
		if(!(model.getEQmode())) params.set("<Lives Left>",fmt.format(model.getAveLL()));
	
		
		if(model.draw()){
			
			
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

	public void run() {
		
		
		// Setup model
		model = new Damage2D(params);
		model.Initialize();
		
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
		
		// Setup output files
		Damage2D.PrintParams(model.getOutdir()+File.separator+"Params.txt",params);	
		model.WriteDataHeaders();
				
		while(model.Avalanche()){
			Job.animate();
			model.TakeDate();
			if (params.sget("Record").equals("On")){
				// prevents taking multiple pics at a single time step
				model.TakePicture(gridS);
				model.TakePicture(gridL);
			}
		}
		
		params.set("<Lives Left>","Finished");
		Job.animate();
		
		return;
	}

}

package chris.old.ofc.apps;

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
import chris.old.ofc.ForceSf;
import chris.old.ofc.QuenchedDamage;

public class ForceSfApp extends Simulation{

	ForceSf model;
	
	Grid gridS = new Grid ("Stress");
	Grid gridL = new Grid ("Lives");
	
	ColorPalette palette1;
	ColorGradient cGradient;
	
	DecimalFormat fmt = new DecimalFormat("0000.00");
	
	double ScMax;
	Boolean pretime, broken;
	
	public static void main(String[] args) {
		new Control(new ForceSfApp(), "OFC Parameters");
	}
	
public void load(Control c) {
		
		params.add("Data Directory",new DirectoryValue("/Users/cserino/Desktop/"));
		params.add("Mode", new ChoiceValue("Earthquake","Damage"));
		params.add("Random Seed",0);
		params.add("Interaction Shape", new ChoiceValue("Circle","Square","Diamond","All Sites"));
		params.add("Interaction Radius (R)",(int)(30));
		params.add("Lattice Size",1<<8);
		params.add("Boundary Condtions", new ChoiceValue("Periodic","Bordered"));
		params.add("Equilibrate Time",(int)0);
		params.add("Number of Lives",2000000);
		params.add("Life Distribution", new ChoiceValue("Constant","Gaussian", "Step Down"));
		params.add("LD Width",0.0);
		params.add("Failure Stress (\u03C3_f)",2.0);
		params.add("\u03C3_f width",(double)(0));
		params.add("d\u03C3_f",0.0001);
		params.add("Residual Stress (\u03C3_r)",1);
		params.add("\u03C3_r width",0.);
		params.add("Dissipation (\u03B1)",new DoubleValue(0.01,0,1));
		params.add("\u03B1 Width", (double)(0));
		params.add("Animation", new ChoiceValue("Off","On"));
		params.addm("Record", new ChoiceValue("Off","On"));
		params.add("\u03C3_f (t)");
		params.add("Last Avalanche Size");
		params.add("<Lives Left>");
		
		c.frameTogether("OFC Model with Damage", gridS, gridL);
		
	}
	
	public void animate() {

		double AS;
		
		broken = ((AS =  model.getAS()) <= 0);
		
		if(pretime){
			params.set("\u03C3_f (t)",model.getTime(-1));
			params.set("Last Avalanche Size",AS);
		}
		else{
			params.set("\u03C3_f (t)",model.getSfofT());
			params.set("Last Avalanche Size",AS);
			if(!(model.getEQmode())) params.set("<Lives Left>",fmt.format(model.getAveLL()));
		}
		
		if(model.draw()){
			
			
			int L = model.getL();
			int N = L*L;
			
			int[] foo = model.getLives();

			double[] copyStress = model.getStress();
			
			for (int jj=0 ; jj < N ; jj++){
				cGradient.getColor(copyStress[jj],-2,model.getSf0());
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
		
		broken = false;
		
		// Setup model
		model = new ForceSf(params);
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
		QuenchedDamage.PrintParams(model.getOutdir()+File.separator+"Params.txt",params);	
		model.WriteDataHeaders();
		
		// Run the simulation
		pretime = true;
		if(model.getTime(-1) != 0){
			params.set("<Lives Left>","Equilibrating");
			while(model.Equilibrate()){
				Job.animate();
			}
		}
		
		pretime = false;	
		model.resetMode(params);	// resets model to Damage mode if specified by params

		if(params.sget("Mode").equals("Earthquake")) params.set("<Lives Left>","Running");
		
		while(model.Avalanche()){
			Job.animate();
			model.TakeDate();
			if (params.sget("Record").equals("On")){
				// prevents taking multiple pics at a single time step
				model.TakePicture(gridS);
				model.TakePicture(gridL);
			}
			
			if(broken) break;
		}
		
		params.set("<Lives Left>","Finished");
		Job.animate();
		
		return;
	}

}

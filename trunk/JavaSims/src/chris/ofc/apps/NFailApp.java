package chris.ofc.apps;


import java.awt.Color;
import java.io.File;
import scikit.graphics.ColorGradient;
import scikit.graphics.ColorPalette;
import scikit.graphics.dim2.Grid;
import scikit.jobs.Control;
import scikit.jobs.Simulation;
import scikit.jobs.params.ChoiceValue;
import scikit.jobs.params.DirectoryValue;
import scikit.jobs.params.DoubleValue;
import chris.ofc.NfailDamage2D;

public class NFailApp extends Simulation {

	Grid grid1 = new Grid ("Stress Lattice");
	Grid grid2 = new Grid ("Failed Sites");

	NfailDamage2D model;
	double ScMax, rgyr;
	
	ColorPalette palette1;
	ColorGradient smooth;
	
	int EEchk = 0;
	

	public static void main(String[] args) {
		new Control(new NFailApp(), "OFC Model");
	}
	
	public void load(Control c) {
		
		params.add("Data Directory",new DirectoryValue("/Users/cserino/Desktop/"));
		params.add("Random Seed",0);
		params.add("Animation", new ChoiceValue("On","Off"));
		params.add("Lattice Size",1<<9);
		params.add("Number of Lives",1);
		params.add("Life Style", new ChoiceValue("Constant","Flat","Gaussian"));
		params.add("Nlives Width",0.1);
		params.add("Boundary Condtions", new ChoiceValue("Periodic","Bordered"));
		params.add("Critical Stress (\u03C3_c)",4.0);
		params.add("\u03C3_c Noise", new ChoiceValue("Off","On"));	
		params.add("\u03C3_c width",Math.sqrt(Math.sqrt(0.4)));
		params.add("Residual Stress (\u03C3_r)",2.0);
		params.add("\u03C3_r Noise", new ChoiceValue("Off","On"));
		params.add("\u03C3_r width",Math.sqrt(Math.sqrt(2)));
		params.add("Interaction Shape", new ChoiceValue("Circle","Square","Diamond","All Sites"));
		params.add("Interaction Radius (R)",(int)(50));
		params.add("Minimum Interaction Radius (r)",0);
		params.add("Dissipation (\u03B1)",new DoubleValue(0.2,0,1));
		params.add("\u03B1 Noise", new ChoiceValue("On","Off"));
		params.add("\u03B1 Width", 0.05);
		params.addm("Record", new ChoiceValue("Off","On"));
		params.add("Number of Resets");
		params.add("Number of Showers");
			
		c.frameTogether("Stress Lattice in an OFC Model with Damage", grid1, grid2);
		
	}
	
	public void animate() {
						
		params.set("Number of Resets",model.time);
		params.set("Number of Showers",model.showernumber);
		
		if (model.ShowGrid){
		
			int[] foo = new int[model.N];
	
				
			for (int i=0 ; i<model.N ; i++){
				smooth.getColor(model.stress[i],-2,ScMax);
				foo[i]=model.alive[i];
			}
				
			grid1.setColors(smooth);
			grid1.registerData(model.L,model.L,model.stress);
			grid2.registerData(model.L, model.L, foo);
				
			if (params.sget("Record").equals("On")){
				model.TakePicture(grid1);
				model.TakePicture(grid2);
				
			}
		
		}

	}

	public void clear() {

		grid1.clear();
		grid2.clear();
	}

	public void run() {
		
		// Initialize Model
		
		model = new NfailDamage2D(params);
		
		String anmt = params.sget("Animation");
		
		if (anmt.equals("On")){
			model.ShowGrid=true;
		}
		else{
			model.ShowGrid=false;
		}
				
		model.Initialize("Flat");
		
		ScMax=model.GetMax(model.Sc);		

		// Set up color scheme
		
		palette1 = new ColorPalette();
		smooth   = new ColorGradient();
		
		int max = model.GetMax(model.alive);

		Color[] Carray = new Color[]{Color.RED,Color.YELLOW,Color.GREEN,Color.BLUE,Color.GRAY};
		
		palette1.setColor(0,Color.WHITE);
		for (int i = 1 ; i <= max ; i++){
			palette1.setColor(i,Carray[i%5]);
		}
		
		grid2.setColors(palette1);
		
		
		// Set up file
		NfailDamage2D.PrintParams(model.outdir+File.separator+"Params.txt", params);	
		model.WriteDataHeader();
		
		
		while(!(model.crack)) {
			
			model.Avalanche();

			model.TakeData();
			
		}
		
	}

}

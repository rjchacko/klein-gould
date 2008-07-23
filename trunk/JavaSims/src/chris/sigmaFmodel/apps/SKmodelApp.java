package chris.sigmaFmodel.apps;

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
import chris.sigmaFmodel.SKdamageModel;

public class SKmodelApp extends Simulation{
	
	SKdamageModel model;
	
	private DecimalFormat fmtD = new DecimalFormat("00000000");

	Grid gridS = new Grid ("Stress");
	Grid gridF = new Grid ("Failures");
	
	ColorPalette palette1;
	ColorGradient cGradient;
	
	double UB;
	
	boolean draw;
		
	public static void main(String[] args) {
		new Control(new SKmodelApp(), "SK Parameters");
	}
	
	public void load(Control c) {
	
		params.add("Data Directory", new DirectoryValue("/home/cserino/Desktop/"));
		params.add("Seed", (int) 0);
		params.add("Interaction Shape", new ChoiceValue("Circle","Square","Diamond","All Sites"));
		params.add("Interaction Radius (R)", (int)(30));
		params.add("Lattice Size (L)", (int) 1<<8);
		params.add("Boundary Condtions", new ChoiceValue("Periodic","Bordered"));
		params.add("Failure Stress (\u03C3_f)", (double) 20.0);
		params.add("\u03C3_f width", (double)(0));
		params.add("Residual Stress (\u03C3_r)", (double) 19.0);
		params.add("\u03C3_r width", (double)0);
		params.add("d(\u03C3_r)", (double) 0.5);
		params.add("d(\u03C3_r) width", (double) 0);
		params.add("Dissipation (\u03B1)",new DoubleValue(0.01,0,1));
		params.add("\u03B1 width", (double)(0));
		params.add("Animation", new ChoiceValue("Off","On"));
		params.addm("Record", new ChoiceValue("Off","On"));
		params.add("Last Avalanche Size");
		params.add("<stress>");
		
		c.frameTogether("S-K Damage Model", gridS, gridF);
		
		return;
	}

	public void run() {

		params.set("<stress>", "Initializing");
		
		model = new SKdamageModel(params);
		
		setupColors();
		setupData();

		while(model.avalanche()){
			
			model.takeData();
			
			if(params.sget("Record").equals("On")){
				model.takePicture(gridS,draw);
				model.takePicture(gridF,draw);
			}
			
			Job.animate();
		}
		
		params.set("<stress>", "Finished");
		
	}

	public void animate() {
		
		params.set("<stress>", model.getSbar());
		params.set("Last Avalanche Size", fmtD.format(model.getAvlnchSize()));
		
		if(!draw) return;
	
		int L = model.getL();
		int N = L*L;

		double[] foo = model.getFailures();
		double[] copyStress = model.getStress();

		for (int jj=0 ; jj < N ; jj++){
			cGradient.getColor(copyStress[jj],0,UB);
		}

		gridS.setColors(cGradient);
		gridS.registerData(L,L,copyStress);
		gridF.registerData(L, L, foo);

	}

	public void clear() {
		
		gridS.clear();
		gridF.clear(); 
	}
	
	private void setupColors(){
		
		draw = (params.sget("Animation").equals("On"));
		UB = params.fget("Failure Stress (\u03C3_f)") + params.fget("\u03C3_f width");
		
		palette1  = new ColorPalette();
		cGradient = new ColorGradient();
	
		Color[] Carray = new Color[]{Color.YELLOW,Color.BLUE,Color.RED,Color.GREEN,Color.GRAY};		
		palette1.setColor(0,Color.BLACK);
		for (int jj = 1 ; jj < model.intN() ; jj++){
			palette1.setColor(jj,Carray[jj%5]);
		}
		palette1.setColor(0,Color.WHITE);
		palette1.setColor(-1,Color.BLACK);
		gridF.setColors(palette1);
		
		return;
	}
	
	private void setupData(){
		
		SKdamageModel.printParams(model.getOutdir()+File.separator+"Params.txt",params);	
		model.writeDataHeaders();
	}
	
}

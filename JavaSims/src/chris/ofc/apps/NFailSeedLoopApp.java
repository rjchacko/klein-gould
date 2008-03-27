package chris.ofc.apps;


import java.io.File;
import java.text.DecimalFormat;
import chris.ofc.NfailDamage2D;
import scikit.graphics.ColorGradient;
import scikit.graphics.ColorPalette;
import scikit.jobs.Control;
import scikit.jobs.Job;
import scikit.jobs.Simulation;
import scikit.jobs.params.ChoiceValue;
import scikit.jobs.params.DirectoryValue;
import scikit.jobs.params.DoubleValue;
import chris.util.PrintUtil;

public class NFailSeedLoopApp extends Simulation {

	NfailDamage2D model; 
	double ScMax, rgyr;
	DecimalFormat fmt = new DecimalFormat("0000000");
	DecimalFormat fmts = new DecimalFormat("0000");
	
	ColorPalette palette1;
	ColorGradient smooth;
	
	int EEchk = 0;
	int SN, cyclevar;
	

	public static void main(String[] args) {
		new Control(new NFailSeedLoopApp(), "OFC Model");
	}
	
	public void load(Control c) {
		
		params.add("Data Directory",new DirectoryValue("/Users/cserino/CurrentSemester/Research/Data/"));
		params.add("Random Seed",0);
		params.add("Number of Sims", 1);
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
		params.add("Interaction Shape", new ChoiceValue("Circle","Square","Diamond"));
		params.add("Interaction Radius (R)",(int)(50));
		params.add("Minimum Interaction Radius (r)",0);
		params.add("Dissipation (\u03B1)",new DoubleValue(0.2,0,1));
		params.add("\u03B1 Noise", new ChoiceValue("On","Off"));
		params.add("\u03B1 Width", 0.05);
		params.add("Number of Resets");
		params.add("Number of Showers");
		params.add("Step Number");
					
	}
	
	public void animate() {
							
		params.set("Number of Resets",model.time);
		params.set("Number of Showers",model.showernumber);
		params.set("Step Number", cyclevar);

	}

	public void clear() {

	}

	public void run() {

		PrintUtil.printlnToFile(params.sget("Data Directory")+File.separator+"Params.txt",params.toString());
		
		for (cyclevar=0; cyclevar < params.iget("Number of Sims") ; cyclevar++){
		
			model = new NfailDamage2D(params);
			
			model.outfile = model.outdir + File.separator+ "Damage" + fmts.format(cyclevar) +".txt";
			model.Initialize(params.sget("Stress Distribution"));
			PrintUtil.printlnToFile(model.outfile,"Time","N_avlnchs","Rgyr","Omega","<FS_stress>","rho_FS");
		
			while(!(model.crack)) {
				
				model.Avalanche();
				model.TakeData();
				
			}
			
			//update seed
			params.set("Random Seed",params.iget("Random Seed")+1);
			
			Job.animate();

		}
		
		System.out.println("All simulations done.");
	
	}
	
	
}

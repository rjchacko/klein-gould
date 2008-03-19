package chris.ofc.apps;


import java.io.File;
import java.text.DecimalFormat;

import scikit.graphics.ColorGradient;
import scikit.graphics.ColorPalette;
import scikit.jobs.Control;
import scikit.jobs.Job;
import scikit.jobs.Simulation;
import scikit.jobs.params.ChoiceValue;
import scikit.jobs.params.DirectoryValue;
import scikit.jobs.params.DoubleValue;
import chris.ofc.NfailDamage2D;
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
		params.add("Looping Parameter", new ChoiceValue("LatticeSize","NumberofLives","CriticalStress","ResidualStress","InteractionRadius","DissipationConstant"));
		params.add("Max Value of LP",10.);
		params.add("LP Step Size",2.);
		params.add("Lattice Size",1<<9);
		params.add("Number of Lives",1);
		params.add("Life Style", new ChoiceValue("Constant","Flat","Gaussian"));
		params.add("Nlives Width",0.1);
		params.add("Boundary Condtions", new ChoiceValue("Periodic","Bordered"));
		params.add("Stress Distribution", new ChoiceValue("Flat","Hammer Blow"));
		params.add("Hammer Size",1);	
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
		
		cyclevar=0;
		
		PrintUtil.printlnToFile(params.sget("Data Directory")+File.separator+"Params.txt",params.toString());

		switch(LoopParams.valueOf(params.sget("Looping Parameter"))) {
		
		case LatticeSize:
			
			for (int lpv = params.iget("Lattice Size") ; lpv <= (int)(params.fget("Max Value of LP")) ; lpv+=(int)(params.fget("LP Step Size"))){

				model = new NfailDamage2D(params);
				
				model.outfile = model.outdir + File.separator+ "Damage" + fmts.format(cyclevar++) +".txt";
									
				model.Initialize(params.sget("Stress Distribution"));
				
				// Set up file
					
				PrintUtil.printlnToFile(model.outfile,"Time","N_avlnchs","Rgyr","Omega","<FS_stress>","rho_FS");
				
				while(!(model.crack)) {
					
					model.Avalanche();
	
					TakeData();
					
				}
				
				Job.animate();
			}
			
			break;
			
		case NumberofLives:
			
			for (int lpv = params.iget("Number of Lives") ; lpv <= (int)(params.fget("Max Value of LP")) ; lpv+=(int)(params.fget("LP Step Size"))){

				model = new NfailDamage2D(params);
				
				model.outfile = model.outdir + File.separator+ "Damage" + fmts.format(cyclevar++) +".txt";
									
				model.Initialize(params.sget("Stress Distribution"));
				
				// Set up file
					
				PrintUtil.printlnToFile(model.outfile,"Time","N_avlnchs","Rgyr","Omega","<FS_stress>","rho_FS");
				
				while(!(model.crack)) {
					
					model.Avalanche();
	
					TakeData();
					
				}
				
				Job.animate();
			}
			
			
			break;
			
		case CriticalStress:
			
			for (double lpv = params.iget("Critical Stress (\u03C3_c)") ; lpv <= params.fget("Max Value of LP"); lpv+=params.fget("LP Step Size")){

				model = new NfailDamage2D(params);
				
				model.outfile = model.outdir + File.separator+ "Damage" + fmts.format(cyclevar++) +".txt";
									
				model.Initialize(params.sget("Stress Distribution"));
				
				// Set up file
					
				PrintUtil.printlnToFile(model.outfile,"Time","N_avlnchs","Rgyr","Omega","<FS_stress>","rho_FS");
				
				while(!(model.crack)) {
					
					model.Avalanche();
	
					TakeData();
					
				}
				
				Job.animate();
			}
			
			
			break;
			
		case Residualtress:
			
			for (double lpv = params.iget("Residual Stress (\u03C3_r)") ; lpv <= params.fget("Max Value of LP") ; lpv+=params.fget("LP Step Size")){

				model = new NfailDamage2D(params);
				
				model.outfile = model.outdir + File.separator+ "Damage" + fmts.format(cyclevar++) +".txt";
									
				model.Initialize(params.sget("Stress Distribution"));
				
				// Set up file
					
				PrintUtil.printlnToFile(model.outfile,"Time","N_avlnchs","Rgyr","Omega","<FS_stress>","rho_FS");
				
				while(!(model.crack)) {
					
					model.Avalanche();
	
					TakeData();
					
				}
				
				Job.animate();
			}
			
			
			break;
			
		case InteractionRadius:
			
			for (int lpv = params.iget("Interaction Radius (R)") ; lpv <= (int)(params.fget("Max Value of LP")) ; lpv+=(int)(params.fget("LP Step Size"))){

				model = new NfailDamage2D(params);
				
				model.outfile = model.outdir + File.separator+ "Damage" + fmts.format(cyclevar++) +".txt";
									
				model.Initialize(params.sget("Stress Distribution"));
				
				// Set up file
					
				PrintUtil.printlnToFile(model.outfile,"Time","N_avlnchs","Rgyr","Omega","<FS_stress>","rho_FS");
				
				while(!(model.crack)) {
					
					model.Avalanche();
	
					TakeData();
					
				}
				
				Job.animate();
			}
			
			
			break;
			
		case DissipationConstant:
			
			for (double lpv = params.iget("Dissipation (\u03B1)") ; lpv <= params.fget("Max Value of LP") ; lpv+=params.fget("LP Step Size")){

				model = new NfailDamage2D(params);
				
				model.outfile = model.outdir + File.separator+ "Damage" + fmts.format(cyclevar++) +".txt";
									
				model.Initialize(params.sget("Stress Distribution"));
				
				// Set up file
					
				PrintUtil.printlnToFile(model.outfile,"Time","N_avlnchs","Rgyr","Omega","<FS_stress>","rho_FS");
				
				while(!(model.crack)) {
					
					model.Avalanche();
	
					TakeData();
					
				}
				
				Job.animate();
			}
			
			
			break;
			
		default:
			System.err.println("Loop Parameter Not Recognized.");
			break;
		}

		System.out.println("Job Finished.");
		
	}
	
	public void TakeData(){
	
		int[] LS = model.LiveSites(); 
		if(LS.length>0){
			rgyr=model.radiusGyration(LS[model.rand.nextInt(LS.length)]);
		}
		else{
			rgyr=0;
		}
		
		PrintUtil.printlnToFile(model.outfile,model.time,model.Nshowers,rgyr,model.EFmetric(),model.GetAve(model.SonFS,model.SonFSindex),model.DaboutFS);
				
		return;
	}
	
	public enum LoopParams {
		LatticeSize, NumberofLives, CriticalStress, Residualtress, InteractionRadius, DissipationConstant;
	}
	
}

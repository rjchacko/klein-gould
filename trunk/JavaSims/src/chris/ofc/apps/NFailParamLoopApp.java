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

public class NFailParamLoopApp extends Simulation {

	NfailDamage2D model; 
	double ScMax, rgyr;
	DecimalFormat fmt = new DecimalFormat("0000000");
	DecimalFormat fmts = new DecimalFormat("0000");
	
	ColorPalette palette1;
	ColorGradient smooth;
	
	int EEchk = 0;
	int SN, cyclevar;
	

	public static void main(String[] args) {
		new Control(new NFailParamLoopApp(), "OFC Model");
	}
	
	public void load(Control c) {
		
		params.add("Data Directory",new DirectoryValue("/Users/cserino/CurrentSemester/Research/Data/"));
		params.add("Random Seed",0);
		params.add("Looping Parameter (LP)", new ChoiceValue("LatticeSize","Nlives","CritStress","ResidStress","InteractionRange","DissipationConstant"));
		params.add("LP Max Value",1<<12);	
		params.add("LP Step Size Value",100);	
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
		
		cyclevar=0;
		
		//PrintUtil.printlnToFile(params.sget("Data Directory")+File.separator+"Params.txt",params.toString());
		model.PrintParams(params.sget("Data Directory")+File.separator+"Params.txt", params);	

		
		switch (LVtype.valueOf(params.sget("Looping Parameter (LP)"))){
		
		case LatticeSize:
			
			for (int lpvar = params.iget("Lattice Size") ; lpvar <= (int)(params.fget("LP Max Value")) ; lpvar+=(int)(params.fget("LP Step Size Value"))){
			
				model = new NfailDamage2D(params);
				
				params.set("Lattice Size",lpvar);
				
				model.outfile = model.outdir + File.separator+ "Damage" + fmts.format(cyclevar++) +".txt";
									
				model.Initialize("Flat");
				
				// Set up file
					
				PrintUtil.printlnToFile(model.outfile,"Lattice Size = ", params.iget("Lattice Size"));
				model.WriteDataHeader(model.outfile);

				
				while(!(model.crack)) {
					
					model.Avalanche();
		
					model.TakeData();
					
				}
				
				Job.animate();
			}
			
			break;
			
		case Nlives:
			
			for (int lpvar = params.iget("Number of Lives") ; lpvar <= (int)(params.fget("LP Max Value")) ; lpvar+=(int)(params.fget("LP Step Size Value"))){
				
				model = new NfailDamage2D(params);
				
				model.outfile = model.outdir + File.separator+ "Damage" + fmts.format(cyclevar++) +".txt";
							
				params.set("Number of Lives",lpvar);
				
				model.Initialize("Flat");
				
				// Set up file
				PrintUtil.printlnToFile(model.outfile,"Number of Lives = ", params.iget("Number of Lives"));
				model.WriteDataHeader(model.outfile);

				while(!(model.crack)) {
					
					model.Avalanche();
		
					model.TakeData();
					
				}
				
				Job.animate();
			}
			
			break;
			
		case CritStress:
			
			for (double lpvar = params.fget("Critical Stress (\u03C3_c)") ; lpvar <= (params.fget("LP Max Value")) ; lpvar+=(params.fget("LP Step Size Value"))){
				
				model = new NfailDamage2D(params);
				
				params.set("Critical Stress (\u03C3_c)",lpvar);
				
				model.outfile = model.outdir + File.separator+ "Damage" + fmts.format(cyclevar++) +".txt";
									
				model.Initialize("Flat");
				
				// Set up file
				
				PrintUtil.printlnToFile(model.outfile,"Critical Stress (\u03C3_c) = ", params.fget("Critical Stress (\u03C3_c)"));
				model.WriteDataHeader(model.outfile);

				while(!(model.crack)) {
					
					model.Avalanche();
		
					model.TakeData();
					
				}
				
				Job.animate();
			}
			
			break;
			
		case ResidStress:
			
			for (double lpvar = params.fget("Residual Stress (\u03C3_r)") ; lpvar <= (params.fget("LP Max Value")) ; lpvar+=(params.fget("LP Step Size Value"))){
				
				model = new NfailDamage2D(params);
				
				model.outfile = model.outdir + File.separator+ "Damage" + fmts.format(cyclevar++) +".txt";
							
				params.set("Residual Stress (\u03C3_r)",lpvar);
				
				model.Initialize("Flat");
				
				// Set up file
				
				PrintUtil.printlnToFile(model.outfile,"Residual Stress (\u03C3_r) = ", params.fget("Residual Stress (\u03C3_r)"));
				model.WriteDataHeader(model.outfile);

				while(!(model.crack)) {
					
					model.Avalanche();
		
					model.TakeData();
					
				}
				
				Job.animate();
			}
			
			break;
			
		case InteractionRange:
			
			for (int lpvar = params.iget("Interaction Radius (R)") ; lpvar <= (int)(params.fget("LP Max Value")) ; lpvar+=(int)(params.fget("LP Step Size Value"))){
				
				model = new NfailDamage2D(params);
				
				model.outfile = model.outdir + File.separator+ "Damage" + fmts.format(cyclevar++) +".txt";
						
				params.set("Interaction Radius (R)",lpvar);
				
				model.Initialize("Flat");
				
				// Set up file
				
				PrintUtil.printlnToFile(model.outfile,"Interaction Radius (R) = ", params.iget("Interaction Radius (R)"));
				model.WriteDataHeader(model.outfile);
				
				while(!(model.crack)) {
					
					model.Avalanche();
		
					model.TakeData();
					
				}
				
				Job.animate();
			}
			
			break;
			
		case DissipationConstant:
			
			for (double lpvar = params.fget("Dissipation (\u03B1)") ; lpvar <= (params.fget("LP Max Value")) ; lpvar+=(params.fget("LP Step Size Value"))){
				
				model = new NfailDamage2D(params);
				
				model.outfile = model.outdir + File.separator+ "Damage" + fmts.format(cyclevar++) +".txt";
									
				params.set("Dissipation (\u03B1)",lpvar);
				
				model.Initialize("Flat");
				
				// Set up file
				
				PrintUtil.printlnToFile(model.outfile,"Dissipation (\u03B1) = ", params.fget("Dissipation (\u03B1)"));
				model.WriteDataHeader(model.outfile);
				
				while(!(model.crack)) {
					
					model.Avalanche();
		
					model.TakeData();
					
				}
				
				Job.animate();
			}
			
			break;
			
		default:
			
			System.err.println("Loop Variable " + params.sget("Looping Parameter (LP)") + " not applicable.");
			
			break;
			

		
		
		}
		
	}
	
	public enum LVtype {
		LatticeSize, Nlives, CritStress, ResidStress, InteractionRange, DissipationConstant;
	}
	
}

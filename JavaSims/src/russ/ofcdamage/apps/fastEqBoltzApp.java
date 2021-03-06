package russ.ofcdamage.apps;

import java.io.File;
import java.text.DecimalFormat;

import russ.ofcdamage2.ofc2Dfast;
import scikit.jobs.Control;
import scikit.jobs.Job;
import scikit.jobs.Simulation;
import scikit.jobs.params.ChoiceValue;
import scikit.jobs.params.DirectoryValue;

public class fastEqBoltzApp extends Simulation{

	private int simt, eqt;
	private ofc2Dfast model;
	private static DecimalFormat fmtI = new DecimalFormat("000");
	
	public static void main(String[] args) {
		new Control(new fastEqBoltzApp(), "OFC Parameters");
	}
	
	public void load(Control c) {
		
		params.add("Data Directory",new DirectoryValue("/Users/cserino/Desktop/"));
		params.add("Data File", "default");
		params.add("Random Seed", (int) 0);
		params.add("Interaction Shape", new ChoiceValue("Circle","Square","Diamond","All Sites"));
		params.add("Interaction Radius (R)", (int) 10);
		params.add("Lattice Size", (int) 256);
		params.add("Boundary Condtions", new ChoiceValue("Periodic","Open"));
		params.add("Equil Time", 100000);
		params.add("Sim Time", 100000);
		params.add("Failure Stress (\u03C3_f)", 2.);
		params.add("\u03C3_f width", 0.);
		params.add("Residual Stress (\u03C3_r)", 1.);
		params.add("\u03C3_r width", 0.05);
		params.add("Dissipation (\u03B1)", 0.05);
		params.add("\u03B1 width", 0.);
		params.add("Status");
		
	}
	
	public void run() {
		
		// Setup model
		params.set("Status", "Intializing");
		Job.animate();
		model = new ofc2Dfast(params);
		model.PrintParams(model.getOutdir()+File.separator+"Params_"+model.getBname()+".txt",params);	
		eqt   = params.iget("Equil Time");
		simt  = params.iget("Sim Time");
		params.set("Status", "Ready");
		Job.animate();
		
		for (int kk = 0 ; kk < 100 ; kk++){
			
			model.setBname(params.sget("Data File")+fmtI.format(kk));
			
			// Equilibrate the system
			for (int jj = 0 ; jj < eqt ; jj++){
				model.evolve(jj,false);
				if(jj%500 == 0){
					params.set("Status", (jj-eqt));
					Job.animate();
				}
			}

			// Simulate the model without damage
			for (int jj = 0 ; jj < simt ; jj++){
				model.evolve(jj,true);
				if(jj%500 == 0){
					params.set("Status", jj*(kk+1));
					Job.animate();
				}
			}

			if((simt-1)%ofc2Dfast.dlength != 0) model.writeData(simt);

			model.clearData();			
		}
		params.set("Status", "Done");
		Job.animate();
		
		return;
	}

	public void animate() {
		
		return;
	}

	public void clear() {
		
		return;
	}

	
}

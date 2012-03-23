package russ.ofcdamage2.apps;

import java.io.File;

import russ.ofcdamage2.ofc2Dfast;
import scikit.jobs.Control;
import scikit.jobs.Job;
import scikit.jobs.Simulation;
import scikit.jobs.params.ChoiceValue;
import scikit.jobs.params.DirectoryValue;

public class fastEqApp extends Simulation{

	private int simt, eqt;
	private ofc2Dfast model;


	public static void main(String[] args) {
		new Control(new fastEqApp(), "OFC Parameters");
	}
	
	public void load(Control c) {
		
		params.add("Data Directory",new DirectoryValue("/Users/cserino/Desktop/"));
		params.add("Data File", "default");
		params.add("Random Seed", (int) 0);
		params.add("Interaction Shape", new ChoiceValue("Circle","Square","Diamond","All Sites"));
		params.add("Interaction Radius (R)", (int) 10);
		params.add("Lattice Size", (int) 256);
		params.add("Boundary Condtions", new ChoiceValue("Open","Periodic"));
		params.add("Equil Time", 100000);
		params.add("Sim Time", 100000);
		params.add("Failure Stress (\u03C3_f)", 2.);
		params.add("\u03C3_f width", 0.);
		params.add("Residual Stress (\u03C3_r)", 1.);
		params.add("\u03C3_r width", 0.025);
		params.add("Dissipation (\u03B1)", 0.025);
		params.add("\u03B1 width", 0.);
		params.add("Status");
		
	}
	
	public void run() {
			
			params.set("Random Seed", params.iget("Random Seed") + 1);
			params.set("Status", "Intializing");
			Job.animate();
			model = new ofc2Dfast(params);
			model.PrintParams(model.getOutdir()+File.separator+"Params_"+model.getBname()+".log",params,this);	
			eqt    = params.iget("Equil Time");
			simt   = params.iget("Sim Time");
			params.set("Status", "Ready");
			Job.animate();

			// Equilibrate the system
			for (int jj = 0 ; jj < eqt ; jj++){
				model.evolve(jj,false);
				if(jj%500 == 0){
					params.set("Status", (jj-eqt));
					Job.animate();
				}
			}

			// Simulate the model in earthquake mode
			for (int jj = 1 ; jj < simt ; jj++){
				model.evolve(jj,true);
				if(jj%500 == 0){
					params.set("Status", jj);
					Job.animate();
				}
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

package russ.ofcdamage.apps;

import java.io.File;

import chris.util.PrintUtil;
import russ.ofcdamage2.ofc2Dfast;
import scikit.jobs.Control;
import scikit.jobs.Job;
import scikit.jobs.params.ChoiceValue;
import scikit.jobs.params.DirectoryValue;
import scikit.jobs.Simulation;

public class startCatalogueApp extends Simulation{

	private ofc2Dfast model;
	
	public static void main(String[] args) {
		new Control(new startCatalogueApp(), "OFC Parameters");
	}
	
	public void load(Control c) {
		
		params.add("Data Directory",new DirectoryValue("/Users/cserino/Documents/Catalogue"));
		params.add("Data File", "default");
		params.add("Random Seed", (int) 0);
		params.add("Interaction Shape", new ChoiceValue("Circle","Square","Diamond","All Sites"));
		params.add("Interaction Radius (R)", (int) 20);
		params.add("Lattice Size", (int) 512);
		params.add("Boundary Condtions", new ChoiceValue("Open","Periodic"));
		params.add("Equil Time", 1000000);
		params.add("Sim Time", 0);
		params.add("Failure Stress (\u03C3_f)", 2.);
		params.add("\u03C3_f width", 0.);
		params.add("Residual Stress (\u03C3_r)", 1.);
		params.add("\u03C3_r width", 0.025);
		params.add("Dissipation (\u03B1)");
		params.add("\u03B1 width", 0.);
		params.add("Status");
	}

	public void run() {

		double[] av = new double[]{0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9};
		int teq  = params.iget("Equil Time");
		for(double alpha : av){
			params.set("Status", "Intializing");
			Job.animate();
			params.set("Dissipation (\u03B1)",alpha);
			model = new ofc2Dfast(params);
			if(alpha == av[0])
				model.PrintParams(model.getOutdir()+File.separator+"Params_"+model.getBname()+".log",params,model);	
			params.set("Status", "Ready");
			Job.animate();
			for (int tt = 0; tt < teq; tt++){
				model.evolve(tt,false);
				if(tt%1000 == 0){
					params.set("Status", (tt-teq));
					Job.animate();
				}
			}
			
			PrintUtil.printVectorToFile(model.getOutdir()+File.separator+model.getBname()+"_Stress_"+100*alpha+".txt",model.getStress());
			params.set("Random Seed", params.iget("Random Seed")+1);
		}



	}
	
	public void animate() {return;}
	public void clear() {return;}
	
}

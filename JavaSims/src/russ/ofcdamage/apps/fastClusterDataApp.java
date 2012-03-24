package russ.ofcdamage.apps;


import java.io.File;
import java.text.DecimalFormat;

import russ.ofcdamage2.damage2Dfast;
import scikit.dataset.Histogram;
import scikit.jobs.Control;
import scikit.jobs.Job;
import scikit.jobs.Simulation;
import scikit.jobs.params.DirectoryValue;
import scikit.jobs.params.Parameters;
import chris.util.PrintUtil;
import chris.util.dummyParamUtil;

public class fastClusterDataApp extends Simulation{

	private int simt, eqt;
	private damage2Dfast model;
	private DecimalFormat pfmt = new DecimalFormat("0000");
	
	private Histogram hes;
	
	public static void main(String[] args) {
		new Control(new fastClusterDataApp(), "OFC Parameters");
	}
	
	public void load(Control c) {
		
		params.add("Data Directory",new DirectoryValue("/Users/cserino/Desktop/"));
		params.add("Data File", "default");
		params.add("Random Seed", (int) 0);
//		params.add("Interaction Shape", new ChoiceValue("Circle","Square","Diamond","All Sites"));
		params.add("Interaction Radius (R)", (int) 10);
		params.add("Lattice Size", (int) 256);
//		params.add("Boundary Condtions", new ChoiceValue("Periodic","Open"));
		params.add("Equil Time", 100000);
		params.add("N_events / Phi", 100000);
		params.add("d(phi)", 0.05);
		params.add("Failure Stress (\u03C3_f)", 2.);
//		params.add("\u03C3_f width", 0.);
		params.add("Residual Stress (\u03C3_r)", 1.);
		params.add("\u03C3_r width", 0.1);
		params.add("Dissipation (\u03B1)", 0.05);
		params.add("\u03D5");
		params.add("Mode");
		params.add("Status");
		
	}
	
	public void run() {
		
		double phin = 1.;
		double dphi = params.fget("d(phi)"); 
		
		// Setup model
		eqt   = params.iget("Equil Time");
		simt  = params.iget("N_events / Phi");
		params.set("Status", "Ready");
		Job.animate();
		
		while (phin > 0){

			params.set("Mode", "Freeze");
			params.set("\u03D5",phin);
			Parameters dparams = dummyParamUtil.ofcParams(params);
			params.set("Status", "Intializing");
			Job.animate();
			hes   = new Histogram(1.);
			model = null;
			model = new damage2Dfast(dparams);
			if(phin == 1){
				model.PrintParams(model.getOutdir()+File.separator+"Params_"+model.getBname()+".log",params);
				PrintUtil.printlnToFile(model.getOutdir()+File.separator+"Params_"+model.getBname()+".log","d(phi) = ",dphi);
			}
			params.set("Status", "Ready");
			Job.animate();
			
			// Equilibrate the system
			for (int jj = 0 ; jj < eqt ; jj++){
				model.evolve(jj,false);
				if(jj%1000 == 0){
					params.set("Status", (jj-eqt));
					Job.animate();
				}
			}

			// Simulate the model without damage
			for (int jj = 0 ; jj < simt ; jj++){
				model.evolve(jj,false);
				hes.accum(model.getGR());
				if(jj%1000 == 0){
					params.set("Status", jj);
					Job.animate();
				}
			}

			printHist(phin);
			phin -= dphi;
			phin = (double)(Math.round(100*phin))/100;
			params.set("Random Seed",params.iget("Random Seed")+1);  
		
			params.set("Status", "Next ");
			Job.animate();
		
		}
		return;
	}

	public void printHist(double phi){
		PrintUtil.printHistToFile(model.getOutdir()+File.separator+model.getBname()+"_"+pfmt.format(100*phi)+".txt",hes);
		hes.clear();
	}
	
	public void animate() {
		
		return;
	}

	public void clear() {
		
		return;
	}

	
}

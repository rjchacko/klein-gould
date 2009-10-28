package chris.ofcdamage.apps;


import java.io.File;
import java.text.DecimalFormat;

import scikit.dataset.Histogram;
import scikit.jobs.Control;
import scikit.jobs.Job;
import scikit.jobs.Simulation;
import scikit.jobs.params.DirectoryValue;
import scikit.jobs.params.Parameters;
import chris.ofcdamage.ofc2Dfast;
import chris.util.PrintUtil;
import chris.util.dummyParamUtil;

public class fastAlphaEffApp extends Simulation{

	private int simt, eqt;
	private ofc2Dfast model;
	private DecimalFormat pfmt = new DecimalFormat("000");
	
	private Histogram hes = new Histogram(1.);
	
	public static void main(String[] args) {
		new Control(new fastAlphaEffApp(), "OFC Parameters");
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
		params.add("Sim Time", 100000);
		params.add("Failure Stress (\u03C3_f)", 2.);
//		params.add("\u03C3_f width", 0.);
		params.add("Residual Stress (\u03C3_r)", 1.);
		params.add("\u03C3_r width", 0.1);
		params.add("Dissipation (\u03B1)", 0.05);
//		params.add("\u03B1 width", 0.);
		params.add("Status");
		
	}
	
	public void run() {
		
		double phin = 1;
		double dphi = 0.1; 
		
		// Setup model
		eqt   = params.iget("Equil Time");
		simt  = params.iget("Sim Time");
		params.set("Status", "Ready");
		Job.animate();
		
		while (phin > 0){
		
			// rescale alpha
			double anow = 1-phin*(1-0.05);
				   anow = (double)(Math.round(10000*anow))/10000;
			params.set("Dissipation (\u03B1)",anow);
			Parameters dparams = dummyParamUtil.ofcParams(params);
			params.set("Status", "Intializing");
			Job.animate();
			model = new ofc2Dfast(dparams);
			if(phin == 1) model.PrintParams(model.getOutdir()+File.separator+"Params_"+model.getBname()+".log",params);	
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
			//if((simt-1)%ofc2Dfast.dlength != 0) model.writeData(simt);

			params.set("Status", "Done");
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

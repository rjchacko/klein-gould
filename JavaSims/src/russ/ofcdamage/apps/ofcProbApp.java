package russ.ofcdamage.apps;

import java.io.File;
import scikit.dataset.Histogram;
import scikit.jobs.Control;
import scikit.jobs.Job;
import scikit.jobs.Simulation;
import scikit.jobs.params.ChoiceValue;
import scikit.jobs.params.DirectoryValue;
import chris.ofcdamage.ofc2Dfast;
import chris.util.PrintUtil;

public class ofcProbApp extends Simulation{

	private int simt, eqt;
	private ofc2Dfast model;
	private Histogram hs, hx;
	
	public static void main(String[] args) {
		new Control(new ofcProbApp(), "OFC Parameters");
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
		params.add("Residual Stress (\u03C3_r)");
		params.add("\u03C3_r width", 0.05);
		params.add("Dissipation (\u03B1)", 0.05);
		params.add("\u03B1 width", 0.);
		params.add("Status");
		
	}
	
	public void run() {
		
		// Setup model
		params.set("Status", "Intializing");
		params.set("Residual Stress (\u03C3_r)", params.fget("\u03C3_r width"));
		Job.animate();
		model = new ofc2Dfast(params);
		model.PrintParams(model.getOutdir()+File.separator+"Params_"+model.getBname()+".txt",params);	
		eqt   = params.iget("Equil Time");
		simt  = params.iget("Sim Time");
		hs    = new Histogram(1e-3);
		hx    = new Histogram(1e-3);
		
		//should be able to get an array of about 1e8 elements (for 16 byte integer / long)
		
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
		
		// Simulate the model without damage
		for (int jj = 0 ; jj < simt ; jj++){
			model.evolve(jj,false);
			takePSdata(hs);
			//model.takeWdata(hx);
			if(jj%500 == 0){
				params.set("Status", jj);
				Job.animate();
			}
		}

		if((simt-1)%ofc2Dfast.dlength != 0) model.writeData(simt);
		printPdata();
				
		params.set("Status", "Done");
		Job.animate();
		
		return;
	}

	private void takePSdata(Histogram h){

		h.accum(model.getStress(0));
		return;
	}
	
	public void animate() {
		
		return;
	}

	public void clear() {
		
		return;
	}

	public void printPdata(){
		
		PrintUtil.printHistToFile(model.getOutdir()+File.separator+model.getBname()+"_Shist"+".txt", hs);
		PrintUtil.printHistToFile(model.getOutdir()+File.separator+model.getBname()+"_Whist"+".txt", hx);
		return;
	}
	
}

package chris.ofcdamage.apps;

import java.io.File;

import scikit.dataset.Histogram;
import scikit.jobs.Control;
import scikit.jobs.Job;
import scikit.jobs.Simulation;
import scikit.jobs.params.ChoiceValue;
import scikit.jobs.params.DirectoryValue;
import chris.ofcdamage.ofc2Dfast;
import chris.util.PrintUtil;

public class fastStrainRateApp extends Simulation{

	private int simt, eqt;
	private ofc2Dfast model;
	public Histogram hGR;

	public static void main(String[] args) {
		new Control(new fastStrainRateApp(), "OFC Parameters");
	}
	
	public void load(Control c) {
		
		params.add("Data Directory",new DirectoryValue("/Users/cserino/Desktop/"));
		params.add("Data File", "default");
		params.add("Random Seed", (int) 0);
		params.add("d\u025B/dt", (double) 2e-5);
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

		double sr;
		model = new ofc2Dfast(params);
		model.PrintParams(model.getOutdir()+File.separator+"Params_"+model.getBname()+".log",params);	
		model = null;
//		for (int pwr = 6 ; pwr > 0 ; pwr--){
//			for (int pf = 1 ; pf < 10 ; pf++){
				
				int pwr = 2;
				int pf  = 1;
				// Setup model
				params.set("Status", "Intializing");
				Job.animate();
				sr     = pf*Math.pow(10, -pwr);
				params.set("d\u025B/dt",sr);
				params.set("Random Seed",params.iget("Random Seed")+1);
				hGR    = new Histogram(1.);
				model  = new ofc2Dfast(params);
				eqt    = params.iget("Equil Time");
				simt   = params.iget("Sim Time");
				params.set("Status", "Ready");
				Job.animate();
				
				// Equilibrate the system
				for (int jj = 0 ; jj < eqt ; jj++){
					model.evolve(jj,false,sr);
					if(jj%500 == 0){
						params.set("Status", (jj-eqt));
						Job.animate();
					}
				}
				
				// Simulate the model without damage
				for (int jj = 0 ; jj < simt ; jj++){
					model.evolve(jj,false,sr);
					hGR.accum(model.getGR());
					if(jj%500 == 0){
						params.set("Status", jj);
						Job.animate();
					}
				}
				
				PrintUtil.printHistToFile(model.getOutdir()+File.separator+"hGR_"+pwr+"_"+pf+".txt", hGR);
				params.set("Status", "Done");
				Job.animate();
				model = null;
//			}
//		}
		
		return;
	}

	public void animate() {
		
		return;
	}

	public void clear() {
		
		return;
	}
}

package chris.ofcdamage.apps;

import java.io.File;
import java.text.DecimalFormat;

import scikit.dataset.Histogram;
import scikit.jobs.Control;
import scikit.jobs.Job;
import scikit.jobs.Simulation;
import scikit.jobs.params.ChoiceValue;
import scikit.jobs.params.DirectoryValue;
import chris.ofcdamage.damage2Dfast;
import chris.util.PrintUtil;

public class grScalingApp2 extends Simulation{

	private int simt, eqt;
	private DecimalFormat pfmt = new DecimalFormat("000");
	private Histogram hist;
	private damage2Dfast model;
	
//	private Grid gridS, gridD;

	public static void main(String[] args) {
		new Control(new grScalingApp2(), "Damage Parameters");
	}
	
	public void load(Control c) {
		
		params.add("Data Directory",new DirectoryValue("/Users/cserino/Desktop/"));
		params.add("Data File", "default");
		params.add("Random Seed", (int) 0);
		params.add("Interaction Shape", new ChoiceValue("Circle","Square","Diamond","All Sites"));
		params.add("Interaction Radius (R)", (int) 10);
		params.add("Lattice Size", (int) 256);
		params.add("Boundary Condtions", new ChoiceValue("Periodic","Open"));
		params.add("Equil Time", 500000);
		params.add("Sim Time", 500000);
		params.add("Number of Lives",(int) 1);
		params.add("NL width", (int) 0);
		params.add("Failure Stress (\u03C3_f)", 2.);
		params.add("\u03C3_f width", 0.);
		params.add("Residual Stress (\u03C3_r)", 1.);
		params.add("\u03C3_r width", 0.1);
		params.add("Dissipation (\u03B1)", 0.05);
		params.add("\u03B1 width", 0.);
		params.add("Mode");
		params.add("\u03D5");
		params.add("Status");
		
//		gridS = new Grid("Stress");
//		gridD = new Grid("Damage");
//		c.frameTogether("Damage Model",gridS,gridD);
		
		return;
	}
	
	public void run() {
		
		double phin = 1;
		simt = params.iget("Sim Time");
		eqt  = params.iget("Equil Time");
		
		while (phin > 0){

			params.set("Mode", "GRscaling");
			params.set("\u03D5",1-phin);
			model  = new damage2Dfast(params);
			hist   = new Histogram(1.);
			
			params.set("Mode", "Equilibrating");
			params.set("\u03D5",phin);
			Job.animate();
			// equilibrate model
			for (int jj = 0 ; jj < eqt ; jj++){
				model.evolveEQ(jj,false);
				//model.evolve(jj,false);
				if(jj%1000 == 0){
					params.set("Status", (jj-eqt));
					Job.animate();
				}
			}

			params.set("Mode", "Simulating");
			params.set("\u03D5",phin);
			Job.animate();
			// simulate model for data
			for (int jj = 0 ; jj < simt ; jj++){
				model.evolveEQ(jj,false);
				//model.evolve(jj,false);
				hist.accum(model.getGR());
				if(jj%1000 == 0){
					params.set("Status", (jj));
					Job.animate();
				}
			}

			writeGRdata(phin);
			phin -= 0.01;
			phin = (double)(Math.round(100*phin))/100;
		}
		params.set("Mode", "Done");
		Job.animate();
		return;
	}

	public void animate() {
		
		return;
	}

	public void clear() {

		return;
	}
	
//	public void animate() {
//		
//		int L = model.getL();
//		
//		params.set("GR Size",model.getGR());
//		gridS.registerData(L, L, model.getStress());
//		gridD.registerData(L, L, model.getDorA());
//
//		return;
//	}
//
//	public void clear() {
//			
//		gridS.clear();
//		gridD.clear();
//		return;
//	}

	private void writeGRdata(double label){
		
		PrintUtil.printHistToFile(model.getOutdir()+File.separator+model.getBname()+"_"+pfmt.format(100*label)+".txt",hist);
		return;
	}

	
}

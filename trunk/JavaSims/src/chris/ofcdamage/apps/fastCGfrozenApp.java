package chris.ofcdamage.apps;

import java.io.File;
import java.text.DecimalFormat;

import scikit.jobs.Control;
import scikit.jobs.Job;
import scikit.jobs.Simulation;
import scikit.jobs.params.ChoiceValue;
import scikit.jobs.params.DirectoryValue;
import chris.ofcdamage.damage2DfastCG;
import chris.ofcdamage.ofc2Dfast;

public class fastCGfrozenApp extends Simulation{

	private int simt, eqt;
	private damage2DfastCG model;
	private DecimalFormat pfmt = new DecimalFormat("000");
	
	public static void main(String[] args) {
		new Control(new fastCGfrozenApp(), "OFC Parameters");
	}
	
	public void load(Control c) {
		
		params.add("Data Directory",new DirectoryValue("/Users/cserino/Desktop/"));
		params.add("Data File", "default");
		params.add("Random Seed", (int) 0);
		params.add("Interaction Shape", new ChoiceValue("Circle","Square","Diamond","All Sites"));
		params.add("Interaction Radius (R)", (int) 10);
		params.add("Lattice Size", (int) 125);
		params.add("L_cg (must be odd)", 5);
		params.add("Boundary Condtions", new ChoiceValue("Periodic","Open"));
		params.add("Number of Lives",(int) 1);
		params.add("NL width", (int) 0);
		params.add("Equil Time", 100000);
		params.add("Sim Time", 100000);
		params.add("t_cg", 100);
		params.add("Failure Stress (\u03C3_f)", 2.);
		params.add("\u03C3_f width", 0.);
		params.add("Residual Stress (\u03C3_r)", 1.);
		params.add("\u03C3_r width", 0.05);
		params.add("Dissipation (\u03B1)", 0.05);
		params.add("\u03B1 width", 0.);
		params.add("Status");
		params.add("\u03D5");
		params.add("Mode");

	}
	
	public void run() {

		double phin = 1;

		eqt  = params.iget("Equil Time");
		simt = params.iget("Sim Time");
		
		// Setup model
		params.set("Status", "Intializing");
		Job.animate();
		
		while (phin > 0){
			
			params.set("Mode", "Freeze");
			params.set("\u03D5",1-phin);
			model  = new damage2DfastCG(params);
			if(phin == 1) model.PrintParams(model.getOutdir()+File.separator+"Params_"+model.getBname()+".log",params,model);	
		
			// Equilibrate the system
			params.set("Mode","Equilibrating");
			params.set("\u03D5",phin);
			Job.animate();
			for (int jj = 0 ; jj < eqt ; jj++){
				model.evolveEQ(jj,false);
				if(jj%500 == 0){
					params.set("Status", (jj-eqt));
					Job.animate();
				}
			}

			// Simulate the model
			params.set("Mode", "Simulating");
			params.set("\u03D5",phin);
			model.setBname(params.sget("Data File")+pfmt.format(100*phin));
			Job.animate();
			for (int jj = 0 ; jj < simt ; jj++){
				model.evolveEQ(jj,true);
				if(jj%500 == 0){
					params.set("Status", jj);
					Job.animate();
				}
			}
			if((simt-1)%ofc2Dfast.dlength != 0) model.writeData(simt);

			//phin -= 0.01;
			phin -= 0.1;
			phin = (double)(Math.round(100*phin))/100;
		}
		
		return;
	}

	public void animate() {
		
		return;
	}

	public void clear() {
		
		return;
	}

	
}

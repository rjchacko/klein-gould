package chris.ofcdamage.apps;

import java.io.File;
import java.text.DecimalFormat;

import scikit.jobs.Control;
import scikit.jobs.Job;
import scikit.jobs.Simulation;
import scikit.jobs.params.ChoiceValue;
import scikit.jobs.params.DirectoryValue;
import chris.ofcdamage.damage2Dfast;
import chris.ofcdamage.damage2DfastCG;

public class fastDamageCgApp extends Simulation{

	private int simt, eqt, Ndead, N;
	private double phis;
	private damage2DfastCG model;
	private DecimalFormat cfmt = new DecimalFormat("00");
	
	public static void main(String[] args) {
		new Control(new fastDamageCgApp(), "OFC Parameters");
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
//		params.add("Number of Lives",(int) 10);
//		params.add("NL width", (int) 8);
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
		params.add("Dead Sites");
		params.add("Mode");

	}
	
	public void run() {

		int dmt  = simt;
		
		// Setup model
		params.set("Status", "Intializing");
		Job.animate();
		model = new damage2DfastCG(params);
		model.PrintParams(model.getOutdir()+File.separator+"Params_"+model.getBname()+".txt",params);	
		eqt   = params.iget("Equil Time");
		simt  = params.iget("Sim Time");
		N     = model.getN();
		Ndead = 0;
		phis  = 0.05;
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

		params.set("Mode","Earthquake");
		// Simulate the model without damage
		for (int jj = 0 ; jj < simt ; jj++){
			model.evolve(jj,true);
			if(jj%500 == 0){
				params.set("Status", jj);
				Job.animate();
			}
		}

		// Simulate the model with damage

		while(Ndead < N){
			
			// simulate with damage until phi >= phis
			params.set("Mode","Damage");
			while(Ndead < phis*N){
				Ndead = model.evolveD(dmt,false);
				if(dmt%500 == 0){
					params.set("Status", (dmt));
				}
				params.set("Dead Sites", Ndead);
				Job.animate();
				dmt++;
				if(Ndead >= N) break;
				System.out.print(Ndead);
				System.out.print("\t");
				System.out.print(phis*N);
			}

			// write and clear data
			if((dmt-1)%damage2Dfast.dlength != 0) model.writeData(dmt);

			if(Ndead >= N) break;

			model.clearData();

			// change save file name
			model.setBname(params.sget("Data File")+"_"+cfmt.format(100*phis));
			//model.PrintParams(model.getOutdir()+File.separator+"Params_"+model.getBname()+".txt",params);

			// re-equilibrate
			params.set("Mode","Equilibrate");
			for (int jj = 0 ; jj < eqt ; jj++){
				model.evolve(jj,false);
				if(jj%500 == 0){
					params.set("Status", (jj-eqt));
					Job.animate();
				}
			}

			// simulate system in EQ mode with phi != 0
			params.set("Mode","Earthquake");
			for (int jj = 0 ; jj < simt ; jj++){
				model.evolve(jj,true);
				if(jj%500 == 0){
					params.set("Status", jj);
					Job.animate();
				}
			}

			phis += 0.05;
			dmt    = simt;				
		}

		if((dmt-1)%damage2Dfast.dlength != 0) model.writeData(dmt);

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

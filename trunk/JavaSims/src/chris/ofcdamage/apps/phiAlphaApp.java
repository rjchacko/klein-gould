package chris.ofcdamage.apps;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;

import scikit.jobs.Control;
import scikit.jobs.Job;
import scikit.jobs.Simulation;
import scikit.jobs.params.ChoiceValue;
import scikit.jobs.params.DirectoryValue;
import chris.ofcdamage.damage2Dfast;

public class phiAlphaApp extends Simulation{

	private int eqt, dmt, Ndead, L, N;
	private damage2Dfast model;

	public static void main(String[] args) {
		new Control(new phiAlphaApp(), "Damage Parameters");
	}
	
	public void load(Control c) {
		
		params.add("Data Directory",new DirectoryValue("/Users/cserino/Desktop"));
		params.add("Data File", "default");
		params.add("Random Seed", (int) 0);
		params.add("Interaction Shape", new ChoiceValue("Circle","Square","Diamond","All Sites"));
		params.add("Interaction Radius (R)", (int) 10);
		params.add("Lattice Size", (int) 256);
		params.add("Boundary Condtions", new ChoiceValue("Periodic","Open"));
		params.add("Equil Time", 100000);
		//params.add("Equil Time", 500000);
		params.add("Number of Lives");
		params.set("Number of Lives",(int) 1);
		params.add("NL width");
		params.set("NL width", (int) 0);
		params.add("Failure Stress (\u03C3_f)", 2.);
		params.add("\u03C3_f width");
		params.set("\u03C3_f width", 0.);
		params.add("Residual Stress (\u03C3_r)", 1.);
		params.add("\u03C3_r width",0.1);
		//params.set("\u03C3_r width", 0.);
		params.add("\u03B1 width");
		params.set("\u03B1 width", 0.);
		params.add("Animate");
		params.set("Animate", "Off");

		params.add("Dissipation (\u03B1)");
		params.add("Mode");
		params.set("Mode","Initializing");
		params.add("Dead Sites");
		return;
	}
	
	public void run() {
		
		int cycle     = 0;
		int cmax      = 50;
		double alphaC = 0.1;	// TEMP 
		double[] phid = new double[cmax];
		double[] htf  = new double[cmax];
		double[] phi  = new double[2];
		
		while (alphaC < 1){
			cycle = 0;
			while (cycle < cmax){

				// Setup model
				params.set("Mode", "Intializing");
				params.set("Dead Sites", "-");
				Job.animate();
				L      = params.iget("Lattice Size");
				N      = L*L;
				params.set("Dissipation (\u03B1)",alphaC);
				model  = new damage2Dfast(params);
				eqt    = params.iget("Equil Time");
				dmt    = 0;
				Ndead  = 0;
				params.set("Mode", "Ready");
				Job.animate();

				// Equilibrate the system
				for (int jj = 0 ; jj < eqt ; jj++){
					model.evolve(jj,false);
					if(jj%500 == 0){
						params.set("Mode", (jj-eqt));
						Job.animate();
					}
				}

				// Simulate the model with damage
				while(Ndead < N){
					Ndead  = model.evolveD(dmt,false);
					phi[1] = phi[0];
					phi[0] = 1.- (double)(Ndead)/((double)(N));;
					if(dmt%500 == 0) params.set("Mode", (dmt));
					params.set("Dead Sites", Ndead);
					Job.animate();
					dmt++;
				}
				params.set("Mode", "Done");
				Job.animate();
				htf[cycle] = dmt;
				if(phi[0] == 0){
					phid[cycle++] = phi[1];
				}
				else{
					phid[cycle++] = phi[0];
				}
				
				params.set("Random Seed",params.iget("Random Seed") + 1);
			}
				// process and write data
			processData(model.getOutdir()+File.separator+model.getBname(), phid, alphaC, htf);
			alphaC = (double)(Math.round(100*(alphaC+0.2)))/100;
		}
		
		return;
	}

	private void processData(String fout, double[] p, double ac, double[] tf){

		
		try{
			File file = new File(fout+"_"+100*ac+"_M"+".txt");
			PrintWriter pw = new PrintWriter(new FileWriter(file, true), true);
			pw.println("phi \t  t \t alpha = "+100*ac);
			for(int jj = 0 ; jj < p.length ; jj++){
				pw.print(p[jj]);
				pw.print("\t");
				pw.println(tf[jj]);
			}
	
		}
		catch (IOException ex){
			ex.printStackTrace();
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

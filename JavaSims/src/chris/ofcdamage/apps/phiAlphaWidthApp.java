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

public class phiAlphaWidthApp extends Simulation{

	private int eqt, dmt, Ndead, L, N;
	private damage2Dfast model;

	public static void main(String[] args) {
		new Control(new phiAlphaWidthApp(), "Damage Parameters");
	}
	
	public void load(Control c) {
		
		params.add("Data Directory",new DirectoryValue("/Users/cserino/Desktop"));
		params.add("Data File", "default");
		params.add("Random Seed", (int) 0);
		params.add("Interaction Shape");
		params.set("Interaction Shape","All Sites");
		params.add("Interaction Radius (R)");
		params.set("Interaction Radius (R)",10);
		params.add("Lattice Size", (int) 128);
		params.add("Boundary Condtions", new ChoiceValue("Periodic","Open"));
		params.add("Equil Time", 10000);
		//params.add("Equil Time", 500000);
		params.add("Number of Lives");
		params.set("Number of Lives",(int) 1);
		params.add("NL width");
		params.set("NL width", (int) 0);
		params.add("Failure Stress (\u03C3_f)", 2.);
		params.add("\u03C3_f width");
		params.set("\u03C3_f width", 0.);
		params.add("Residual Stress (\u03C3_r)", 1.);
		params.add("\u03C3_r width",0.);
		//params.set("\u03C3_r width", 0.);
		params.add("\u03B1 width");
		params.set("\u03B1 width", 0.);
		params.add("Animate");
		params.set("Animate", "Off");

		params.add("Dissipation (\u03B1)", 0.2);
		params.add("Status (out of 1000)");
		params.add("time");
		return;
	}
	
	public void run() {
		
		int cycle      = -1;
		double[] phi   = new double[2];
		int cn         = 100;
		double[] store = new double[cn];

		while (cycle++ < (cn-1)){

			// Setup model
			params.set("Status (out of 1000)", "Intializing");
			Job.animate();
			L      = params.iget("Lattice Size");
			N      = L*L;
			//params.set("Dissipation (\u03B1)",alphaC);
			model  = new damage2Dfast(params);
			//model.setBname(bn+"_"+fmt.format(100*alphaC)+"_"+fmt.format(cycle));
			if(cycle == 0){
				model.PrintParams(model.getOutdir()+File.separator+"Params_"+model.getBname()+".txt",params);	
				setupOF(model.getOutdir()+File.separator+model.getBname()+"_M"+".txt");
			}
			eqt    = params.iget("Equil Time");
			dmt    = 0;
			Ndead  = 0;
			params.set("Status (out of 1000)", cycle);
			Job.animate();

			// Equilibrate the system
			for (int jj = 0 ; jj < eqt ; jj++){
				model.evolve(jj,false);
				if(jj%1000 == 0){
					params.set("time", (jj-eqt));
					Job.animate();
				}
			}

			// Simulate the model with damage
			while(Ndead < N){
				Ndead  = model.evolveD(dmt,false);
				phi[1] = phi[0];
				phi[0] = 1.- (double)(Ndead)/((double)(N));;
				if(dmt%1000 == 0) params.set("time", (dmt));
				dmt++;
			}
			params.set("time", "Intializing");
			Job.animate();
			// process and write data
			if(phi[0] == 0){
				store[cycle] = phi[1];
			}
			else{
				store[cycle] = phi[0];
			}
			params.set("Random Seed",params.iget("Random Seed") + 1);
		}
		
		processData(model.getOutdir()+File.separator+model.getBname()+"_M"+".txt", store);

		params.set("time", "Done");
		Job.animate();
		
		return;
	}

	private void processData(String fout, double[] p){

		try{
			File file = new File(fout);
			PrintWriter pw = new PrintWriter(new FileWriter(file, true), true);
			for(int jj = 0 ; jj <  p.length ; jj++){
				pw.println(p[jj]);
			}
			pw.close();
		}
		catch (IOException ex){
			ex.printStackTrace();
		}

		return;
	}
	
	private void setupOF(String fout){
		try{
			File file = new File(fout);
			PrintWriter pw = new PrintWriter(new FileWriter(file, true), true);
			pw.println("phi(t_fail)");
			pw.close();
		}
		catch (IOException ex){
			ex.printStackTrace();
		}
	}
	
	
	public void animate() {

		return;
	}

	public void clear() {

		return;
	}
	
}

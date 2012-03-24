package russ.ofcdamage.apps;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;

import russ.ofcdamage2.damage2Dfast;
import scikit.jobs.Control;
import scikit.jobs.Job;
import scikit.jobs.Simulation;
import scikit.jobs.params.ChoiceValue;
import scikit.jobs.params.DirectoryValue;
import chris.util.dummyParamUtil;

public class fastT2Fapp extends Simulation{

	private int eqt;
	private damage2Dfast model;
	
	public static void main(String[] args) {
		new Control(new fastT2Fapp(), "Damage Parameters");
	}
	
	public void load(Control c) {
		
		params.add("Data Directory",new DirectoryValue("/Users/cserino/Desktop/"));
		params.add("Data File", "default");
		params.add("Random Seed", (int) 0);
		params.add("Interaction Radius (R)", (int) 10);
		params.add("Lattice Size", (int) 256);
		params.add("Boundary Condtions", new ChoiceValue("Open","Periodic"));
		params.add("Equil Time", 100000);
		params.add("Number of Lives",(int) 10);
		params.add("NL width", (int) 0);
		params.add("Failure Stress (\u03C3_f)", 2.);
		params.add("Residual Stress (\u03C3_r)", 1.);
		params.add("\u03C3_r width", 0.05);
		params.add("Dissipation (\u03B1)", 0.05);
		params.add("Status");
		params.add("Mode");
		params.add("Iteration");

		return;
	}
	
	public void run() {
		
		eqt = params.iget("Equil Time");
		int cycle = 1;
		
		while (true){
			
			// Setup model
			params.set("Status", "Intializing");
			params.set("Mode", "-");
			model = new damage2Dfast(dummyParamUtil.ofcParams(params));
			if(cycle == 1) model.PrintParams(model.getOutdir()+File.separator+"Params_"+model.getBname()+".txt",params);	

			// Equilibrate the system
			params.set("Mode", "Equilibrating");
			params.set("Iteration", cycle);	
			Job.animate();
			for (int jj = 0 ; jj < eqt ; jj++){
				model.evolve(jj,false);
				if(jj%1000 == 0){
					params.set("Status", (jj-eqt));
					Job.animate();
				}
			}

			// Simulate the model to failure
			params.set("Mode", "Damage");
			Job.animate();
			int Nd = 0;
			int N  = model.getN();
			int t  = 0;
			while(Nd < N){
				t++;
				Nd = model.evolveD(t,false);
				if(t%1000 == 0){
					params.set("Status", Nd);
					Job.animate();
				}
			}
			
			// record data and update random seed
			writeT(t);
			cycle++;
			params.set("Random Seed",params.iget("Random Seed")+1);
		}
	}

	public void animate() {

		return;
	}

	public void clear() {

		return;
	}
	
	private void writeT(int t){
		try{
			File file = new File(model.getOutdir()+File.separator+model.getBname()+"_T2F.txt");
			PrintWriter pw = new PrintWriter(new FileWriter(file, true), true);
			pw.println(t);
			pw.close();
		}
		catch (IOException ex){
			ex.printStackTrace();
		}

		return;
	}
}

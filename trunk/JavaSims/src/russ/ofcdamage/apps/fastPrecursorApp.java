package russ.ofcdamage.apps;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.text.DecimalFormat;

import scikit.jobs.Control;
import scikit.jobs.Job;
import scikit.jobs.Simulation;
import scikit.jobs.params.ChoiceValue;
import scikit.jobs.params.DirectoryValue;
import chris.ofcdamage.damage2Dfast;
import chris.util.dummyParamUtil;

public class fastPrecursorApp extends Simulation{

	private int simt, eqt, Ndead, L, N, gr[];
	private damage2Dfast model;
	private DecimalFormat fmt = new DecimalFormat("00");

	public static void main(String[] args) {
		new Control(new fastPrecursorApp(), "Damage Parameters");
	}
	
	public void load(Control c) {
		
		params.add("Data Directory",new DirectoryValue("/Users/cserino/Desktop/"));
		params.add("Data File", "default");
		params.add("Random Seed", (int) 0);
		params.add("Interaction Shape", new ChoiceValue("Circle","Square","Diamond","All Sites"));
		params.add("Interaction Radius (R)", (int) 10);
		params.add("Lattice Size", (int) 256);
		params.add("Boundary Condtions", new ChoiceValue("Open","Periodic"));
		params.add("Equil Time", 100000);
		params.add("Sim Time", 100000);
		params.add("Number of Lives",(int) 5);
		params.add("Failure Stress (\u03C3_f)", 2.);
		params.add("Residual Stress (\u03C3_r)", 1.);
		params.add("\u03C3_r width", 0.025);
		params.add("Dissipation (\u03B1)", 0.05);
		params.add("Mode");
		params.add("Dead Sites");

		return;
	}
	
	public void run() {
				
		// Setup model
		double anow = 0.05;
		double da   = 0.05;
		
		while(anow < 1){
		
			for (int cycle = 0 ; cycle < 10 ; cycle++){

				params.set("Mode", "Intializing");
				params.set("Dead Sites", "-");
				params.set("Dissipation (\u03B1)", anow);
				Job.animate();
				L      = params.iget("Lattice Size");
				N      = L*L;
				model  = new damage2Dfast(dummyParamUtil.ofcParams(params));
				if(cycle == 0 && anow == 0.05)
					model.PrintParams(model.getOutdir()+File.separator+"Params_"+model.getBname()+".log",params);	
				eqt    = params.iget("Equil Time");
				simt   = params.iget("Sim Time");
				Ndead  = 0;
				gr     = new int[simt+N*params.iget("Number of Lives")];
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

				// Simulate the model without damage
				for (int jj = 0 ; jj < simt ; jj++){
					model.evolve(jj,false);
					gr[jj] = model.getGR();
					if(jj%500 == 0){
						params.set("Mode", jj);
					}
					Job.animate();
				}

				// Simulate the model with damage
				int t = simt;
				while(Ndead < N){
					Ndead = model.evolveD(t,true);
					gr[t] = model.getGR();
					if((t++)%500 == 0){
						params.set("Mode", t);
					}
					params.set("Dead Sites", Ndead);
					Job.animate();
				}

				params.set("Mode", "Writing Data File");
				Job.animate();
				printData(cycle, anow, gr);
				model = null;
				params.set("Random Seed", params.iget("Random Seed")+1);
			}

			params.set("Mode", "Next Alpha");
			Job.animate();
			anow += da;
		}
		params.set("Mode", "Done");
		Job.animate();
		return;
	}

	public void printData(int c, double anow, int[] h){
		try{
			File file = new File(model.getOutdir()+File.separator+model.getBname()+"_"+fmt.format(100*anow)+"_"+fmt.format(c)+".txt");
			PrintWriter pw = new PrintWriter(new FileWriter(file, true), true);
			for (int jj = 0 ; jj < h.length; jj++){
				pw.println(h[jj]);
			}
		} catch (IOException ex){
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

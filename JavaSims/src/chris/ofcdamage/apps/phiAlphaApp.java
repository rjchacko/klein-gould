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
		
		params.add("Data Directory",new DirectoryValue("/Users/cserino/Desktop/Testing"));
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
		params.add("\u03C3_r width");
		params.set("\u03C3_r width", 0.);
		params.add("\u03B1 width");
		params.set("\u03B1 width", 0.);
		params.add("Animate");
		params.set("Animate", "Off");

		params.add("Dissipation (\u03B1)");
		params.add("Status");
		params.add("Dead Sites");
		return;
	}
	
	public void run() {
		
		int cycle     = 0;
		double alphaC = 0;
		double[] phid = new double[5];
		double[] phi  = new double[2];
		
		while (alphaC < 1){
			alphaC = (double)(Math.round(100*(alphaC+0.05)))/100;
			cycle = 0;
			while (cycle < 5){

				// Setup model
				params.set("Status", "Intializing");
				params.set("Dead Sites", "-");
				Job.animate();
				L      = params.iget("Lattice Size");
				N      = L*L;
				params.set("Dissipation (\u03B1)",alphaC);
				model  = new damage2Dfast(params);
				//model.setBname(bn+"_"+fmt.format(100*alphaC)+"_"+fmt.format(cycle));
				if(cycle == 0 && alphaC == 0.05){
					model.PrintParams(model.getOutdir()+File.separator+"Params_"+model.getBname()+".txt",params);	
					setupOF(model.getOutdir()+File.separator+model.getBname()+"_M"+".txt");
				}
				eqt    = params.iget("Equil Time");
				dmt    = 0;
				Ndead  = 0;
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

				// Simulate the model with damage
				while(Ndead < N){
					Ndead  = model.evolveD(dmt,false);
					phi[1] = phi[0];
					phi[0] = 1.- (double)(Ndead)/((double)(N));;
					if(dmt%500 == 0) params.set("Status", (dmt));
					params.set("Dead Sites", Ndead);
					Job.animate();
					dmt++;
				}
				params.set("Status", "Done");
				Job.animate();
				if(phi[0] == 0){
					phid[cycle++] = phi[1];
				}
				else{
					phid[cycle++] = phi[0];
				}
				params.set("Random Seed",params.iget("Random Seed") + 1);
			}
				// process and write data
			processData(model.getOutdir()+File.separator+model.getBname()+"_M"+".txt", phid, alphaC);
		}
		
		return;
	}

	private void processData(String fout, double[] p, double ac){
		double ave = 0;
		double var = 0;
		
		for (int jj = 0 ; jj < 5 ; jj++){
			ave += p[jj];
		}
		ave = ave/5;
		for (int jj = 0 ; jj < 5 ; jj++){
			var += (ave-p[jj])*(ave-p[jj]);
		}
		var = Math.sqrt(var/5);
		
		try{
			File file = new File(fout);
			PrintWriter pw = new PrintWriter(new FileWriter(file, true), true);
			pw.print(ac);
			pw.print("\t");
			pw.print(ave);
			pw.print("\t");
			pw.print(var);
			pw.print("\t");
			pw.print(p[0]);
			pw.print("\t");
			pw.print(p[1]);
			pw.print("\t");
			pw.print(p[2]);
			pw.print("\t");
			pw.print(p[3]);
			pw.print("\t");
			pw.print(p[4]);
			pw.println();
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
			pw.print("alpha");
			pw.print("\t");
			pw.print("<phi>");
			pw.print("\t");
			pw.print("var{phi}");
			pw.print("\t");
			pw.print("phi_1");
			pw.print("\t");
			pw.print("phi_2");
			pw.print("\t");
			pw.print("phi_3");
			pw.print("\t");
			pw.print("phi_4");
			pw.print("\t");
			pw.print("phi_5");
			pw.println();
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

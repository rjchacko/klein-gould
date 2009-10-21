package chris.tests;

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

public class testFrozen extends Simulation{

	private damage2Dfast model;

	public static void main(String[] args) {
		new Control(new testFrozen(), "TESTING");
	}
	
	
	@Override
	public void animate() {
		// TODO Auto-generated method stub
		
	}

	@Override
	public void clear() {
		// TODO Auto-generated method stub
		
	}

	public void load(Control c) {

		params.add("Data Directory",new DirectoryValue("/Users/cserino/Desktop/"));
		params.add("Data File", "default");
		params.add("Random Seed", (int) 2344);
		params.add("Interaction Shape", new ChoiceValue("Circle","Square","Diamond","All Sites"));
		params.add("Interaction Radius (R)", (int) 1);
		params.add("Lattice Size", (int) 5);
		params.add("L_cg (must be odd)", 5);
		params.add("Boundary Condtions", new ChoiceValue("Periodic","Open"));
		params.add("Number of Lives",(int) 2);
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

		double phin = 0.75;

		
		// Setup model
		params.set("Status", "Intializing");
		Job.animate();
		
		params.set("Status", "Regular");
		Job.animate();
		params.set("Mode", "foobar");
		params.set("\u03D5",1-phin);
		model  = new damage2Dfast(params);	
		printRelevantVars();
		
		params.set("Status", "Frozen");
		Job.animate();
		params.set("Mode", "Freeze");
		params.set("\u03D5",1-phin);
		model  = new damage2Dfast(params);	
		printRelevantVars();

		
		params.set("Status", "Done");
		Job.animate();
		return;
	}
	
	public void printRelevantVars(){
		
		try{
			File file = new File("/Users/cserino/Desktop/checkFrozen.txt");
			PrintWriter pw = new PrintWriter(new FileWriter(file, true), true);
			pw.print("Site");
			pw.print("\t");
			pw.print("Live Nbs");
			pw.print("\t");
			pw.print("Lives");
			pw.print("\t");
			pw.print("Sr");
			pw.print("\t");
			pw.print("Sf");
			pw.print("\t");
			pw.print("Stress");
			pw.print("\t");
			pw.print("Alive");
			pw.println();
			int N = params.iget("Lattice Size");
			N     = N*N;
			for(int jj = 0 ; jj < N ; jj++){
				pw.print(jj);
				pw.print("\t");
				pw.print(model.getLN(jj));
				pw.print("\t");
				pw.print(model.getNL(jj));
				pw.print("\t");
				pw.print(model.getSr(jj));
				pw.print("\t");
				pw.print(model.getSf(jj));
				pw.print("\t");
				pw.print(model.getStress(jj));
				pw.print("\t");
				pw.print(model.isAlive(jj));
				pw.println();
			}
			pw.println("-------------------------------------------------------------------");
			
			
			
			pw.close();
		}
		catch (IOException ex){
			ex.printStackTrace();
		}
		
		return;
	}

}

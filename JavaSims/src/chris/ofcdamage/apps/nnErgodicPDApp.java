package chris.ofcdamage.apps;

import java.io.File;

import scikit.jobs.Control;
import scikit.jobs.Job;
import scikit.jobs.Simulation;
import scikit.jobs.params.DirectoryValue;
import chris.foo.ofc.Damage2D;
import chris.ofcdamage.damage;

public class nnErgodicPDApp extends Simulation {

	private damage model;
	private int pt, ptmax;
	
	public static void main(String[] args) {
		new Control(new nnErgodicPDApp(), "OFC Parameters");
	}
	
	
	public void load(Control c) {
		params.add("Data Directory",new DirectoryValue("/Users/cserino/Research/Data2/ErgPhaseDiagram/NearestNeighbor/"));
		params.add("Random Seed",0);
		params.add("Interaction Shape");
		params.add("Interaction Radius (R)");
		params.add("Lattice Size");
		params.add("Boundary Condtions");
		params.add("Intitial Stess");
		params.add("Min Lives");	
		params.add("Max Lives");
		params.add("Failure Stress (\u03C3_f)");
		params.add("\u03C3_f width");
		params.add("Residual Stress (\u03C3_r)");
		params.add("\u03C3_r width",0.5);
		params.add("Dissipation (\u03B1)");
		params.add("\u03B1 width");
		params.add("Equil Time");
		params.add("Trend Time");
		params.add("Animation");
		params.add("Record");
		params.add("Number of Plate Updates");
		params.add("N_dead");

		params.set("Interaction Radius (R)",(int)(1));
		params.set("Interaction Shape","Diamond");
		params.set("Boundary Condtions","Periodic");
		params.set("Lattice Size",1<<8);
		params.set("Intitial Stess","Random");
		params.set("Min Lives", 20);	
		params.set("Max Lives", 20);
		params.set("Equil Time", (int)(1e6));	// 1e6
		params.set("Trend Time", (int)(1e7));	// 1e7
		params.set("Animation", "Off");
		params.set("Record", "Off");
		params.set("Failure Stress (\u03C3_f)",2.0);
		params.set("\u03C3_f width",(double)(0));
		params.set("Residual Stress (\u03C3_r)",1.25);
		params.set("Dissipation (\u03B1)",0.01);
		params.set("\u03B1 width",0);
	}

	public void run() {
		
		// Setup model
		params.set("N_dead","Initializing");
		model = new damage(params);
		
		// Setup output files
		configOF();
				
		// Equilibrate
		equil();
		// Simulate w/o Damage for Data
		ideal();
		
		params.set("N_dead","Done");
		Job.animate();
		
		return;
	}
	
	private void configOF(){
		
		Damage2D.PrintParams(model.getOutdir()+File.separator+"Params.txt",params);	
		model.writeDataHeaders();
		
		return;
	}
	
	private void equil(){
	
		// Equilibrate
		params.set("N_dead", "Equilibrating");
		pt      = 0;
		ptmax   = params.iget("Equil Time");
		
		model.setEquil(true);
		model.runClocks(false);
		while(pt < ptmax){
			model.equilibrate();
			pt++;
			if(pt%100==0) Job.animate();
		}
		
		return;
	}
	
	private void ideal(){
		
		// Simulate w/o Damage for Data
		params.set("N_dead","Simulating w/o damage");
		model.setEquil(true);
		model.runClocks(true);		
		int pt0 = 0;
		ptmax = params.iget("Trend Time"); 
		pt    = params.iget("Trend Time");
		while(pt0 < ptmax){
			model.avalanche();
			model.takeData();
			pt++;
			pt0++;
			if(pt0%100==0) Job.animate();
		}
		
		return;
	}
	
	public void animate() {
		params.set("Number of Plate Updates", pt - ptmax);
	}

	public void clear() {
			
	}
	


}

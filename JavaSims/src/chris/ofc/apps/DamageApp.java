package chris.ofc.apps;

import scikit.graphics.dim2.Grid;
import scikit.jobs.Control;
import scikit.jobs.Simulation;
import scikit.jobs.params.ChoiceValue;
import scikit.jobs.params.DirectoryValue;
import scikit.jobs.params.DoubleValue;

public class DamageApp extends Simulation{

	Grid gridS = new Grid ("Stress");
	Grid gridL = new Grid ("Lives");

	public static void main(String[] args) {
		new Control(new DamageApp(), "OFC Parameters");
	}
	
public void load(Control c) {
		
		params.add("Data Directory",new DirectoryValue("/Users/cserino/Desktop/"));
		params.add("Random Seed",0);
		params.add("Interaction Shape", new ChoiceValue("Circle","Square","Diamond","All Sites"));
		params.add("Interaction Radius (R)",(int)(50));
		params.add("Lattice Size",1<<9);
		params.add("Boundary Condtions", new ChoiceValue("Periodic","Bordered"));
		params.add("Equilibrate Time",(int)100000);
		params.add("Number of Lives",6);
		params.add("Life Distribution", new ChoiceValue("Constant","Gaussian", "Step Down"));
		params.add("LD Width",0.1);
		params.add("Critical Stress (\u03C3_c)",4.0);
		params.add("\u03C3_c width",(double)(0));
		params.add("Residual Stress (\u03C3_r)",2.0);
		params.add("\u03C3_r width",(double)(0));
		params.add("Dissipation (\u03B1)",new DoubleValue(0.2,0,1));
		params.add("\u03B1 Width", 0.05);
		params.add("Animation", new ChoiceValue("Off","On"));
		params.addm("Record", new ChoiceValue("Off","On"));
		params.add("Number of Plate Updates");
		
		c.frameTogether("OFC Model with Damage", gridS, gridL);
		
	}
	
	public void animate() {
		
	}

	public void clear() {
		
	}

	public void run() {
		
	}

}

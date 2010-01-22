package chris.MD.TwoD.apps;

import static scikit.util.Utilities.asList;

import java.text.DecimalFormat;

import scikit.graphics.dim2.Scene2D;
import scikit.jobs.Control;
import scikit.jobs.Job;
import scikit.jobs.Simulation;
import chris.MD.TwoD.IdealGas;

public class IdealGasApp extends Simulation{

	public Scene2D canvas = new Scene2D("Particles");
	public IdealGas model;
	static DecimalFormat tf = new DecimalFormat("######.00");
	
	public static void main(String[] args) {
		new Control(new IdealGasApp(), "Ideal Gas");
	}
	


	public void animate() {

		params.set("t", tf.format(model.gettime()));
		canvas.setDrawables(asList(model.boundaryDw(), model.particlesDw()));
	}
	
	public void clear() {
		
		canvas.clear();
	}


	public void load(Control c) {
		
		params.add("seed",0);
		params.add("Lx",100);
		params.add("Ly",100);
		params.add("R", 5);
		params.add("M",1);
		params.add("N", 2);
		params.add("T", 5);
		params.add("dt", 1.);
		params.add("d\u03C4",1.);
		params.add("t");
		
		c.frame(canvas);		
	}

	public void run() {
		
		params.add("Initial Conditions","Kinetic Energy");
		params.add("Boundary Conditions","Closed");

		model = new IdealGas(params);
		
		while(true){
			model.step();
			Job.animate();
		}
	}

	
	
}

package chris.MD.TwoD.apps;

import static scikit.util.Utilities.asList;

import java.text.DecimalFormat;

import scikit.graphics.dim2.Scene2D;
import scikit.jobs.Control;
import scikit.jobs.Job;
import scikit.jobs.Simulation;
import scikit.jobs.params.ChoiceValue;
import chris.MD.TwoD.LennardJones;

public class LJapp extends Simulation{

	public Scene2D canvas = new Scene2D("Particles");
	public LennardJones model;
	static DecimalFormat tf = new DecimalFormat("######.00");
	static DecimalFormat Ef = new DecimalFormat("0.###E0");
	private double now, then, tau, T;

	
	public static void main(String[] args) {
//		new Control(new LJapp(), "Lennard-Jones System").getJob().throttleAnimation(true);
		new Control(new LJapp(), "Lennard-Jones System");

	}
	
	public void animate() {
		
		double[] eng = new double[2];
		if(T != params.fget("T")){
			T = params.fget("T");
			model.changeT(T);
		}
		tau   = params.fget("d\u03C4");
		params.set("t", tf.format(now));
		model.getEsys(eng);
		params.set("T_m", Ef.format(eng[0]/model.N));
		params.set("E",Ef.format(eng[0]+eng[1]));

		canvas.setDrawables(asList(model.boundaryDw(), model.particlesDw()));
		return;
	}

	public void clear() {

		canvas.clear();
		return;
	}

	public void load(Control c) {
		


		params.add("seed",0); 
		params.add("L",14.);
		params.add("Boundary Conditions", new ChoiceValue("Periodic","Closed"));
		params.add("Initial Conditions", new ChoiceValue("Melt","Read In", "Copy", "Viscous"));
		params.add("ODE Solver", new ChoiceValue("Velocity Verlet","First Order Euler"));
		params.add("N",225);
		params.add("M",1);
		params.add("R",0.25); 
		params.add("dt",5e-3);
		params.addm("d\u03C4",1.);
		params.addm("T",1e-6);
		params.add("t");
		params.add("E");
		params.add("T_m");

		c.frame(canvas);				
	}

	public void run() {

		model = new LennardJones(params);
		tau   = params.fget("d\u03C4");
		T     = params.fget("T");
		then  = -2*tau;
		Job.animate();

		while(true){
			now = model.stepWT();
			if(now - then >= tau){
				then = now;
				Job.animate();
			}
		}
		
	}

}

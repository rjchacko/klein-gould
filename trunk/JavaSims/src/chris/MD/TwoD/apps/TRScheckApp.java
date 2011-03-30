package chris.MD.TwoD.apps;

import static scikit.util.Utilities.asList;

import java.text.DecimalFormat;

import scikit.graphics.dim2.Scene2D;
import scikit.jobs.Control;
import scikit.jobs.Job;
import scikit.jobs.Simulation;
import scikit.jobs.params.ChoiceValue;
import chris.MD.TwoD.TRS;
import chris.util.PrintUtil;

public class TRScheckApp extends Simulation{

	public Scene2D canvas = new Scene2D("Particles");
	public TRS model;
	static DecimalFormat tf = new DecimalFormat("######.00");
	static DecimalFormat Ef = new DecimalFormat("0.###E0");
	private double now, then, tau, T;
	private boolean quench;
	
	public static void main(String[] args) {
		new Control(new TRScheckApp(), "T-R-S System").getJob().throttleAnimation(true);
//		new Control(new TRScheckApp(), "T-R-S System");

	}
	
	public void animate() {
		
		/*
		 * a better thermometer may be 
		 * to keep an average over the past
		 * 10 (or so) time steps?
		 * 
		 */
		
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
		
		if(flags.contains("Quench")){
			flags.clear();
			if(quench)
				return;
			quench = true;
			model.quenchA2(0.013);
		}
		
		
//		PrintUtil.printlnToFile("/Users/cserino/Desktop/EofT.txt",now,eng[0],eng[1]);
		return;
	}

	public void clear() {
		
		canvas.clear();
		return;
	}

	public void load(Control c) {
		
		params.add("seed",0); 
		double lx = 20*1.272596758935476;
		params.add("Lx",lx); // for the perfect hex in debug mode 
		double ly = lx*Math.sqrt(3)/2.;
		params.add("Ly",ly); //  "
//		params.add("Lx",30.); // for the perfect hex in debug mode 
//		params.add("Ly",30.); //  "
		params.add("Boundary Conditions", new ChoiceValue("Periodic","Closed"));
		params.add("Initial Conditions", new ChoiceValue("Melt","Read In", "Copy", "Viscous","Debug"));
		params.add("ODE Solver", new ChoiceValue("Velocity Verlet","First Order Euler"));
		params.add("N",400);
		params.add("M",1);
		params.add("R",0.5); 
		params.add("dt",1e-2);
		params.addm("d\u03C4",0.25);
		params.addm("T",0.);
		params.add("t");
		params.add("E");
		params.add("T_m");
		flags.add("Quench");


		c.frame(canvas);	
	}

	public void run() {
		model = new TRS(params);
		PrintUtil.printlnToFile("/Users/cserino/Desktop/Params.txt",params.toString(this.getClass().getName()));
		tau   = params.fget("d\u03C4");
		T     = params.fget("T");
		then  = -2*tau;
		Job.animate();
		quench = false;

		
//		model.quenchA2(0.013);  // quench to A2 = 0.013 < Ac2 (square is stable)
		while(true){
			now = model.stepWT();
			if(now - then >= tau){
				then = now;
				Job.animate();
			}
		}
		
	}

}

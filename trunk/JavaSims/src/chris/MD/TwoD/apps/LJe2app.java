package chris.MD.TwoD.apps;

import java.text.DecimalFormat;
import scikit.graphics.dim2.Scene2D;
import scikit.jobs.Control;
import scikit.jobs.Job;
import scikit.jobs.Simulation;
import scikit.jobs.params.ChoiceValue;
import chris.MD.TwoD.LennardJones;
import chris.util.DirUtil;
import chris.util.PrintUtil;

public class LJe2app extends Simulation{

	public Scene2D canvas = new Scene2D("Particles");
	public LennardJones model;
	static DecimalFormat tf = new DecimalFormat("######.00");
	static DecimalFormat Ef = new DecimalFormat("0.###E0");
	private double tmax, E[];
	
	public static void main(String[] args) {
		new Control(new LJe2app(), "Lennard-Jones System");

	}
	
	public void animate() {

		params.set("t", tf.format(model.gettime()));
		params.set("T_m", Ef.format(model.getTemp()));
		params.set("E",Ef.format(model.getEsys()));
		return;
	}

	public void clear() {

		canvas.clear();
		return;
	}

	public void load(Control c) {

		params.add("seed",0); 
		params.add("Lx",50);
		params.add("Ly",50); 
		params.add("Boundary Conditions", new ChoiceValue("Periodic","Closed"));
		params.add("Initial Conditions", new ChoiceValue("Melt", "Viscous", "Copy", "Read In"));
		params.add("ODE Solver", new ChoiceValue("Velocity Verlet","First Order Euler"));
		params.add("N",16);
		params.add("M",1);
		params.add("R",1.25); 
		params.add("dt");
		params.addm("T",1e-4);
		params.add("t");
		params.add("E");
		params.add("T_m");
	}

	public void run() {
		
		// generated using log spacing of 100 points between 10^-1 and 10^-4 
		double dt[] = new double[]{0.100000000 ,0.074989421 ,0.056234133 ,0.042169650 ,0.031622777 ,0.023713737 ,0.017782794 ,0.013335214 ,0.010000000 ,0.007498942 ,0.005623413 ,0.004216965 ,0.003162278 ,0.002371374 ,0.001778279 ,0.001333521 ,0.001000000 ,0.000749894 ,0.000562341 ,0.000421697 ,0.000316228 ,0.000237137 ,0.000177828 ,0.000133352 ,0.000100000};
		int dtI, countt;
		
		tmax = 1e5;
		dtI  = 0;
		DirUtil.MkDir("/Users/cserino/Desktop/Efluc");
		PrintUtil.printlnToFile("/Users/cserino/Desktop/Efluc/Params.txt",params.toString());
		while (dtI < dt.length){
			E = new double[(int)(tmax)+2];
			params.set("dt",dt[dtI]);
			model = new LennardJones(params);
			countt  = 0;
			Job.animate();
			E[0] = model.getEsys();
			while(model.gettime() < tmax){
				E[(int)(model.gettime())+1] = model.stepWE();
				if(Double.isNaN(E[(int)(model.gettime())+1]))
					break;
				if(countt++ % 100 == 0)
					Job.animate();
			}
			PrintUtil.printlnToFile("/Users/cserino/Desktop/Efluc/E"+dtI+".txt","dt = ", dt[dtI++]);
			PrintUtil.printVectorToFile("/Users/cserino/Desktop/Efluc/E"+dtI+".txt",E,(int)(model.gettime())+1);
		}
		
	}

}

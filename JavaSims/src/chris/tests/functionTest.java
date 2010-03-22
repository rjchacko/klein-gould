package chris.tests;

import java.awt.Color;

import scikit.dataset.PointSet;
import scikit.graphics.dim2.Plot;
import scikit.jobs.Control;
import scikit.jobs.Job;
import scikit.jobs.Simulation;
import scikit.jobs.params.ChoiceValue;
import chris.MD.TwoD.LennardJones;
import chris.util.vector2d;

public class functionTest extends Simulation{

	private Plot Vplot  = new Plot("Potential");
	private Plot Fplot1 = new Plot("Force Num. Der.");
	private Plot Fplot2 = new Plot("Force Hard Coded");
	private PointSet V, F1, F2;
	
	// for L-J potential
	private LennardJones lj;

	
	public static void main(String[] args) {
		new Control(new functionTest(), "Test Potentials and Forces for MD");
	};
	
	public void animate() {
		
		int N = 10000;
		double rmin = 1.999;
		double rmax = 2.001;
		double[] r   = new double[N+1];
		double[] phi = new double[N+1];
		double[] Frc = new double[N+1];
		double[] dph = new double[N+1];

		
		for(int jj = 0 ; jj < N+1 ; jj++){
			r[jj]   =jj*(rmax-rmin)*params.fget("\u03C3")/N + rmin;
			phi[jj] = lj.potential(new vector2d(r[jj],0));
			Frc[jj] = lj.force(new vector2d(r[jj],0)).x;
		}
		for(int jj = 1 ; jj < N ; jj++){
			dph[jj] = -(phi[jj+1]-phi[jj-1])/(r[jj+1]-r[jj-1]);
		}
		dph[0] = dph[1];
		dph[N] = dph[N-1];
		
		//gEnergy = new PointSet(crd,model.gLandscape());
		V  = new PointSet(r,phi);
		F1 = new PointSet(r,dph);
		F2 = new PointSet(r,Frc);

		Vplot.registerLines("Potential", V, Color.BLACK);
		Fplot1.registerLines("Force", F1, Color.BLACK);
		Fplot2.registerLines("Force", F2, Color.BLACK);

		return;
	}

	public void clear() {
		
	}

	public void load(Control c) {
		
		// for L-J potential
		params.add("seed",0); 
		params.add("Lx",50);
		params.add("Ly",50); 
		params.add("Boundary Conditions", new ChoiceValue("Periodic","Closed"));
		params.add("Initial Conditions", new ChoiceValue("Melt", "Viscous", "Copy", "Read In"));
		params.add("ODE Solver", new ChoiceValue("Velocity Verlet","First Order Euler"));
		params.add("N",100);
		params.add("M",1);
		params.add("R",1.); 
		params.add("\u03C3",1.);
		params.add("dt",1e-2);
		params.addm("d\u03C4",1.);
		params.add("t");
		params.add("E");
		params.add("T");

		c.frameTogether("MD Potential and Force", Vplot, Fplot1,Fplot2);
	}

	public void run() {

		lj = new LennardJones(params);
		
		while(true)
			Job.animate();
		
	}

}

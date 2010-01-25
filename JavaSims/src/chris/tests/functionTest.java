package chris.tests;

import java.awt.Color;

import scikit.dataset.PointSet;
import scikit.graphics.dim2.Plot;
import scikit.jobs.Control;
import scikit.jobs.Job;
import scikit.jobs.Simulation;
import chris.MD.TwoD.LennardJones;

public class functionTest extends Simulation{

	private Plot Vplot = new Plot("Potential");
	private Plot Fplot = new Plot("Force");
	private PointSet V, F;
	
	// for L-J potential
	private LennardJones lj;

	
	public static void main(String[] args) {
		new Control(new functionTest(), "Test Potentials and Forces for MD");
	};
	
	public void animate() {
		
		int N = 10000;
		double[] r   = new double[N];
		double[] phi = new double[N];
		double[] Frc = new double[N];
		
		for(int jj = 0 ; jj < N ; jj++){
			r[jj]   = 1.6+ jj*0.8*params.fget("\u03C3")/N ;
			phi[jj] = lj.potential(r[jj]);
			Frc[jj] = lj.force(r[jj],0).x;
		}
		
		//gEnergy = new PointSet(crd,model.gLandscape());
		V = new PointSet(r,phi);
		F = new PointSet(r,Frc);
		Vplot.registerLines("Potential", V, Color.BLACK);
		Fplot.registerLines("Force", F, Color.BLACK);

		return;
	}

	public void clear() {
		
	}

	public void load(Control c) {
		
		// for L-J potential
		params.add("\u03C3",1.);
		params.add("dt",1);
		params.add("N",2);
		params.add("Lx",2);
		params.add("Ly",2); 
		params.add("seed",0); 
		params.add("Boundary Conditions","Periodic");
		params.add("Initial Conditions","Kinetic Energy");
		params.add("M",1);
		params.add("R",1); 

		c.frameTogether("MD Potential and Force", Vplot, Fplot);
	}

	public void run() {

		lj = new LennardJones(params);

		while(true)
			Job.animate();
		
	}

}

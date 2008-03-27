package rachele.ising.dim2.apps;

import static scikit.util.Utilities.format;

import java.awt.Color;
//import static scikit.util.Utilities.frame;
import scikit.dataset.PointSet;
import scikit.graphics.dim2.Plot;
import scikit.jobs.Control;
import scikit.jobs.Job;
import scikit.jobs.Simulation;
import scikit.jobs.params.ChoiceValue;
import rachele.ising.dim1.FieldIsing1D;

public class StripeClumpLinearTheory extends Simulation {
	Plot clumps = new Plot("Clumps");
    FieldIsing1D ising;
    // 2D matrix has n x n dimensions
    // must have n = 2m+1 where m is positive integer
    public int matrixDim = 5;
    public static final double k0 = 4.4934092;
    public int kyInt = 1; //ky = k0 * kyInt;
    public double [] a = new double [matrixDim];
    public double [] v = new double [matrixDim-2];
	
    public void load(Control c) {
    	c.frame(clumps);
    	params.addm("Model", new ChoiceValue("A", "B"));
    	params.addm("Noise", new ChoiceValue("On", "Off"));
    	params.addm("T", 0.85);
    	params.addm("J", -1.0);
    	params.addm("H", 0.04);
    	params.addm("R", 100000);
    	params.add("L/R", 32.0);
    	params.add("R/dx", 4.0);
    	params.add("kR bin-width", 0.1);
    	params.add("Random seed", 0);
    	params.add("Density", -.4);
    	params.add("dt", 0.1);
    	params.add("Time Allocation");
    	params.add("max Write Time", 30.0);
    	params.add("Time Count");
    	params.add("Time");
    	params.add("DENSITY");
    	params.add("Lp");
    	params.add("F");

    	flags.add("Calculate");
    }
	
	@Override
	public void animate() {
		params.set("Time", format(ising.t));
		clumps.registerLines("Field", new PointSet(0, ising.dx, ising.phi), Color.BLACK);
		if(flags.contains("Calculate")) calculate();
			
	}

	public void clear() {
	}

	public void run() {
		ising = new FieldIsing1D(params);
		ising.simulate();
		Job.animate();
		//run the 1D antiferromagnetic ising model until satisfied	
	}

	private void calculate(){
		//find a_n coefficients of matrix
		findCoefficients();
		//double maxEigenvalue = findMaxEigenvalue();
			
	}
	
	private void findCoefficients(){
		for(int i = 0; i <= matrixDim; i ++){
			double sum = 0;
			for(int point = 0; point < ising.Lp; point ++){
				sum += Math.cos(i*k0*point)/(1-Math.pow(ising.phi[point],2));
			}
			double ave = sum/(double)ising.Lp;
			a[i] = ave;
		}
		for(int i = 0; i <= matrixDim - 2; i++ ){
			double V = (Math.sin(i*k0)/(i*k0))*(Math.sin(kyInt*k0)/(kyInt*k0));
			
			v[i] = -V-ising.T*a[0];
		}
	}
	
//	private double findMaxEigenvalue(){
//		double [] eignevalue = new double [matrixDim];
//		//the matrix looks like:
//		
//		//						...
//		//		v2		a1		a2		a3		a4
//		//		a1		v1		a1		a2		a3
//		//...	a2		a1		v0		a1		a2		...
//		//		a3		a2		a1		v1		a1
//		//		a4		a3		a2		a1		v2
//		//					...		
//		double maxEigenvalue = 
//		return 
//	}
}

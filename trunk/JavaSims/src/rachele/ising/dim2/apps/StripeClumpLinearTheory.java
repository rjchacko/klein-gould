package rachele.ising.dim2.apps;

import static scikit.util.Utilities.format;

import java.awt.Color;
//import java.io.DataInputStream;
//import java.io.EOFException;
//import java.io.File;
//import java.io.FileInputStream;
//import java.io.FileNotFoundException;
//import static scikit.util.Utilities.frame;
import scikit.dataset.PointSet;
import scikit.graphics.dim2.Plot;
import scikit.jobs.Control;
import scikit.jobs.Job;
import scikit.jobs.Simulation;
import scikit.jobs.params.ChoiceValue;
import scikit.numerics.Jama.EigenvalueDecomposition;
import scikit.numerics.Jama.Matrix;
import rachele.ising.dim1.FieldIsing1D;

public class StripeClumpLinearTheory extends Simulation {
	Plot clumps = new Plot("Clumps");
    FieldIsing1D ising;
    // 2D matrix has n x n dimensions
    // must have n = 2m+1 where m is positive integer
    public int matrixDim = 11;
    public static final double k0 = 4.4934092;
    public int kyInt = 1; //ky = k0 * kyInt;
    public double [] a = new double [matrixDim];
    public double [] v = new double [matrixDim-2];
	public EigenvalueDecomposition Eig;

    
	public static void main(String[] args) {
		new Control(new StripeClumpLinearTheory(), "Stripe -> Clump");
	}
    
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
		readInConfiguration();
		findCoefficients();
		constructMatrix();
		//double maxEigenvalue = findMaxEigenvalue();
			
	}
	
	public void readInConfiguration(){
//		try{
//			File myFile = new File("../../../research/javaData/configs1d/config");
//			DataInputStream dis = new DataInputStream(new FileInputStream(myFile));
//			int timeIndex, spaceIndex;
//			double phiValue;
//			try{
//				while(true){
//					timeIndex = dis.readInt();
//					dis.readChar();       // throws out the tab
//					spaceIndex =dis.readInt();
//					dis.readChar();       // throws out the tab
//					phiValue = dis.readDouble();
//					dis.readChar();
//					//phi[timeIndex][spaceIndex] = phiValue;
//				}
//			} catch (EOFException e) {
//			}
//
//		} catch (FileNotFoundException e) {
//			System.err.println("FileStreamsTest: " + e);
//		} catch (Exception ex) {
//			ex.printStackTrace();
//		}
	}
	
	private void findCoefficients(){
		for(int i = 0; i < matrixDim; i ++){
			double sum = 0;
			for(int point = 0; point < ising.Lp; point ++){
				sum += Math.cos(i*k0*point)/(1-Math.pow(ising.phi[point],2));
			}
			double ave = sum/(double)ising.Lp;
			a[i] = ave;
		}
		for(int i = 0; i < matrixDim - 2; i++ ){
			double V = (Math.sin(i*k0)/(i*k0))*(Math.sin(kyInt*k0)/(kyInt*k0));
			
			v[i] = -V-ising.T*a[0];
		}
	}
	private void constructMatrix(){
		double [][] matrix = new double [matrixDim][matrixDim];
		for(int i = 0; i < matrixDim; i++){
			int vIndex = Math.abs(i-((matrixDim-1)/2));
			//System.out.println("i = " + i + " vIndex = " + vIndex);
			matrix[i][i] = v[vIndex];
			for (int j = 0; j < i; j++)
				matrix[j][i] = a[i-j];
		}
		//for (int i = 0; i < matrixDim; i ++)
			//System.out.println(matrix[i][j]);
		Matrix A = new Matrix(matrix);	
		Eig = A.eig();
		//the matrix looks like:
		
		//						...
		//		v2		a1		a2		a3		a4
		//		a1		v1		a1		a2		a3
		//...	a2		a1		v0		a1		a2		...
		//		a3		a2		a1		v1		a1
		//		a4		a3		a2		a1		v2
		//					...	
	}
//	private double findMaxEigenvalue(){
//		double [] eignevalue = new double [matrixDim];
	
//		double maxEigenvalue = 
//		return 
//	}
}

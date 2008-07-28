package rachele.ising.dim2.apps;
//import static java.lang.Math.PI;
//import static java.lang.Math.sin;

import java.awt.Color;
import scikit.dataset.PointSet;
import scikit.graphics.dim2.Plot;
import scikit.jobs.Control;
import scikit.jobs.Simulation;
import scikit.jobs.params.ChoiceValue;
import scikit.jobs.params.FileValue;
import scikit.numerics.Jama.EigenvalueDecomposition;
import scikit.numerics.Jama.Matrix;
//import scikit.numerics.fft.FFT1D;
import scikit.numerics.fft.managed.ComplexDoubleFFT;
import scikit.numerics.fft.managed.ComplexDoubleFFT_Mixed;
import scikit.util.*;
import rachele.ising.dim2.FindCoefficients;

public class StripeClumpLinearTheory extends Simulation{
	Plot sf = new Plot ("Structure Fuction");
	Plot sfPhi = new Plot ("Phi0 SF Plot");
	Plot plot = new Plot ("Input Plot");
	Plot phiPlot = new Plot("Phi Plot");
    // 2D matrix has n x n dimensions
    // must have n = 2m+1 where m is positive integer
	public int matrixDim, vDim, Lp;
	public double T;
	public double [] phi, v, phiSF;
	public int kyInt = 1;
	public double dx, L, R, kChunk;
	public EigenvalueDecomposition Eig;
	public double[] fftScratch;
	public FindCoefficients coeff;
    
	public static void main(String[] args) {
		new Control(new StripeClumpLinearTheory(), "Stripe -> Clump");
	}
	
	public void load(Control c) {
		params.add("Input File", new FileValue("/home/erdomi/data/lraim/configs1d/phi0"));
		params.addm("Center On", new ChoiceValue("low frequencies", "high frequencies"));
		params.addm("Matrix Dim", 15);
		params.addm("Temperature", 0.04);
		params.addm("Lp", 128);
		params.addm("R", 2000000);
		params.addm("L/R", 3.0);
		c.frameTogether("f(x) vs phi(x)", plot, sf, phiPlot, sfPhi);
	}
	
	public void calculate(){
		T = params.fget("Temperature");
		matrixDim = params.iget("Matrix Dim");
		Lp = params.iget("Lp");
		R = params.fget("R");
		L = params.fget("L/R")*R;
		dx = L/Lp;
		kChunk = 2*Math.PI*R/L;
		phi = new double [Lp];
		phiSF = new double [Lp];
	
		vDim = 1+(matrixDim-1)/2;
		v = new double [vDim];
		String inputFile = params.sget("Input File");
		coeff = new FindCoefficients(Lp, inputFile);
		findVelements();
		constructMatrix();
		double maxEigenvalue = findMaxEigenvalue();
		System.out.println("matrix dim = " + matrixDim + " final slope = " +2*maxEigenvalue);		
	}

	private void findVelements(){
		for(int i = 0; i < vDim; i++ ){
			double Vx = (i == 0) ? 1.0 :  Math.sin(i*kChunk)/(i*kChunk);
			double Vy = (i == 0) ? 1.0 :  Math.sin(kyInt*kChunk)/(kyInt*kChunk);
			double V = Vx*Vy;
			coeff.findCoefficientsFromFile();
			v[i] = -V-T*coeff.a[0];
		}		
	}
	
	public void constructMatrix(){
		double [][] matrix = new double [matrixDim][matrixDim];
		for(int i = 0; i < matrixDim; i++){
			int vIndex = Math.abs(i-((matrixDim-1)/2));
			//System.out.println("i = " + i + " vIndex = " + vIndex);
			matrix[i][i] = v[vIndex];
			for (int j = 0; j < i; j++)
				matrix[j][i] = matrix [i][j] = -T*coeff.a[i-j];
		}
//		for (int i = 0; i < matrixDim; i ++){
//			for (int j = 0; j < matrixDim; j ++){
//				System.out.println(i + " " + j + " " + matrix[i][j]);	
//			}
//		}
		printMatrix(matrix);
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
		// where vi = v[i] = -V(k)-T a0
		// and ai = -T a[i]	
	}
	
	private void printMatrix(double[][] matrix) {
		System.out.println("Here comes the matrix");
		for (int i = 0; i < matrixDim; i ++){
			for(int j = 0; j < matrixDim; j++){
				System.out.print(matrix[i][j] + "     ");
			}
			System.out.println(" ");
		}
	}

	private double findMaxEigenvalue(){
		double [] eigenvalue = new double [matrixDim];
		eigenvalue = Eig.getRealEigenvalues();
		for (int i = 0; i < matrixDim; i ++)
			System.out.println("eigenvalue " + i + " = " + eigenvalue[i]);
		double maxEigenvalue = DoubleArray.max(eigenvalue);
		Matrix V = Eig.getV();
		Matrix D = Eig.getD();
		double[] finalDvector = D.transpose().getArray()[matrixDim - 1];
		double[] finalEigenvector = V.transpose().getArray()[matrixDim - 1];
		for (int i = 0; i < matrixDim; i ++){
			System.out.println(i + " " + finalEigenvector[i]);
		}
		System.out.println("eigenvector ? " + finalDvector[matrixDim - 1]);
		return maxEigenvalue;
	}


	
	public PointSet getSF(){
		return new PointSet(0, 1, coeff.a);
	}
	
	public PointSet getInput(){
		return new PointSet(0, 1, coeff.inputFunction);
	}

	public PointSet getPhiSlice(){
		return new PointSet(0, 1, coeff.phiSlice);
	}
	
	public PointSet getPhiSF(){
		ComplexDoubleFFT fft = new ComplexDoubleFFT_Mixed(Lp);
		double [] fftScratch = new double [2*Lp];
		for (int i = 0; i < Lp; i++) {
			fftScratch[2*i] = coeff.phiSlice[i];
			fftScratch[2*i+1] = 0;
		}
		fft.transform(fftScratch);
		for (int i = 0; i < Lp; i++) 
			phiSF[i] = fftScratch[2*i];	
		//		FFT1D fft = new FFT1D(Lp);
//		fft.setLength(L);
//		fft.transform(phi, new FFT1D.MapFn() {
//			public void apply(double k, double re, double im) {
//				double kR = k*R;
//				int kRLatticePoint = (int)Math.abs(kR/kChunk);
//				phiSF[kRLatticePoint] = re;
//				//System.out.println("point " + kRLatticePoint + " " + a[kRLatticePoint]);
//			}
//		});		
		return new PointSet(0, 1, phiSF);
	}
	
	public void animate() {
		phiPlot.registerLines("Slice", getPhiSlice(), Color.BLUE);
		sf.registerLines("SF", getSF(), Color.BLACK);
		sfPhi.registerLines("Phi SF",getPhiSF(), Color.RED);
		plot.registerLines("input", getInput(), Color.BLACK);
	}

	public void clear() {
	}

	public void run() {
		calculate();
	}
	
}

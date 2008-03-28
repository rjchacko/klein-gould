package rachele.ising.dim2.apps;
import scikit.graphics.dim2.Plot;
import scikit.jobs.Control;
import scikit.jobs.Simulation;
import scikit.jobs.params.ChoiceValue;
import scikit.numerics.Jama.EigenvalueDecomposition;
import scikit.numerics.Jama.Matrix;
import scikit.util.*;
import rachele.util.FileUtil;

public class StripeClumpLinearTheory extends Simulation{
	Plot clumps = new Plot("Clumps");
    // 2D matrix has n x n dimensions
    // must have n = 2m+1 where m is positive integer
	public int matrixDim;
	public int vDim;
	//input parameters:
	public double T = 0.04;
	public int configSize;
	public double [] phi;
		
	//public static final double k0 = 4.4934092;
	public static final double k0 = 4.188790204;
	public int kyInt = 1; //ky = k0 * kyInt;
    public double [] a; 
    public double [] v; 
	public EigenvalueDecomposition Eig;
    
	public static void main(String[] args) {
		new Control(new StripeClumpLinearTheory(), "Stripe -> Clump");
	}
	
	public void load(Control c) {
		params.addm("Center On", new ChoiceValue("low frequencies", "high frequencies"));
		params.addm("Matrix Dim", 5);
		params.addm("Lp", 128);
		params.addm("R", 2000000);
	}
	
	public void calculate(){
		matrixDim = params.iget("Matrix Dim");
		configSize = params.iget("Lp");
		double kStep = 2*Math.PI*params.fget("R")/configSize;
		System.out.println("kstep = " + kStep);
		phi = new double [configSize];
		a = new double [matrixDim];
		vDim = 1+(matrixDim-1)/2;
		v = new double [vDim];
		//find a_n coefficients of matrix
		String fileName = "../../../research/javaData/configs1d/config";
		phi = FileUtil.readConfigFromFile(fileName, configSize);
//		for (int i = 0; i < configSize; i ++)
//			System.out.println(i + " " + phi[i]);
		findCoefficients(phi, configSize);
		constructMatrix();
		double maxEigenvalue = findMaxEigenvalue();
		System.out.println(matrixDim + "   " +2*maxEigenvalue);		
	}
	
	public void findCoefficients(double [] A, int size){
		for(int i = 0; i < matrixDim; i ++){
			double sum = 0;
			for(int point = 0; point < size; point ++){
				sum += Math.cos(i*k0*point)/(1-Math.pow(A[point],2));
			}
			double ave = sum/(double)size;
			a[i] = ave;
			//System.out.println("a" + i + " = " + a[i]);			
		}
		for(int i = 0; i < vDim; i++ ){
			double Vx;
			if (i ==0)
				Vx = 1.0;
			else
				Vx = Math.sin(i*k0)/(i*k0);
			double Vy;
			if (kyInt == 0)
				Vy = 1.0;
			else 
				Vy = Math.sin(kyInt*k0)/(kyInt*k0);
			double V = Vx*Vy;
			v[i] = -V-T*a[0];
			//System.out.println("v" + i + " = " + v[i]);
		}
	}
	public void constructMatrix(){
		double [][] matrix = new double [matrixDim][matrixDim];
		for(int i = 0; i < matrixDim; i++){
			int vIndex = Math.abs(i-((matrixDim-1)/2));
			//System.out.println("i = " + i + " vIndex = " + vIndex);
			matrix[i][i] = v[vIndex];
			for (int j = 0; j < i; j++)
				matrix[j][i] = matrix [i][j] = -T*a[i-j];
		}
//		for (int i = 0; i < matrixDim; i ++){
//			for (int j = 0; j < matrixDim; j ++){
//				System.out.println(i + " " + j + " " + matrix[i][j]);	
//			}
//		}
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
	private double findMaxEigenvalue(){
		double [] eigenvalue = new double [matrixDim];
		eigenvalue = Eig.getRealEigenvalues();
		for (int i = 0; i < matrixDim; i ++)
			System.out.println("eigenvalue " + i + " = " + eigenvalue[i]);
		double maxEigenvalue = DoubleArray.max(eigenvalue);
		Matrix V = Eig.getV();
		//Matrix D = Eig.getD();
		for (int i = 0; i < matrixDim; i ++){
			for (int j = 0; j < matrixDim; ){
				
			}
		}
		double[] eigenvector0 = V.transpose().getArray()[0];
		double[] eigenvector1 = V.transpose().getArray()[1];
		System.out.println(eigenvector0[0] + " " + eigenvector0[1]);
		System.out.println(eigenvector1[0] + " " + eigenvector1[1]);
		return maxEigenvalue;
	}

	public void animate() {
	}

	public void clear() {
	}

	public void run() {
		calculate();
	}
	
}

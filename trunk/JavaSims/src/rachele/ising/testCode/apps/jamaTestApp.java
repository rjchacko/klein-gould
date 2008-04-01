package rachele.ising.testCode.apps;

import scikit.numerics.Jama.*;

public class jamaTestApp{

    public static void main(String[] argv)
    {     
        double [][] M6 = new double[2][2];
	
	M6[0][0] = 1.0;
	M6[1][0] = 1.0;
    M6[0][1] = 1.0;
	M6[1][1] = 1.0;
	Matrix A = new Matrix(M6);	
	EigenvalueDecomposition Eig = A.eig();
	
	
	Matrix V = Eig.getV();
	double[] eigenvector0 = V.transpose().getArray()[0];
	double[] eigenvector1 = V.transpose().getArray()[1];
	System.out.println(eigenvector0[0] + " " + eigenvector0[1]);
	System.out.println(eigenvector1[0] + " " + eigenvector1[1]);
	
	Matrix D = Eig.getD();
	double[] eigenvalue0 = D.transpose().getArray()[0];
	double[] eigenvalue1 = D.transpose().getArray()[1];
	System.out.println(eigenvalue0[0] + " " + eigenvalue0[1]);
	System.out.println(eigenvalue1[0] + " " + eigenvalue1[1]);

	double[] e = Eig.getRealEigenvalues();
	System.out.println(e[0] + " " + e[1]);
    }
}

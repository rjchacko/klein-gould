package rachele.ising.testCode.apps;

import scikit.numerics.Jama.*;
import rachele.util.MathTools;

public class jamaTestApp{

    public static void main(String[] argv){   
    int L = 1012;
   	double [][] M = new double[L][L];
	
   	
   	for(int x = 0; x < L; x++){
   		for(int y = 0; y <= x; y++){
   			M[x][y] = M[y][x] = x*y;
   		}
   	}

//   	System.out.println("Here comes the matrix:");
//   	for(int x = 0; x < L; x++){
//   		for(int y = 0; y < L; y++){
//   			System.out.print("" + M[x][y] + " ");
//   		}
//   		System.out.println(" ");
//   	}
   	
	Matrix A = new Matrix(M);	
	EigenvalueDecomposition Eig = A.eig();
	
	
	Matrix V = Eig.getV();
	double [][] VV = V.transpose().getArray();//take transpose to get eigenvectors as row

//	System.out.println(VV[0][0] + " " + VV[0][1]+ " " + VV[0][2]);
//	System.out.println(VV[1][0] + " " + VV[1][1]+ " " + VV[1][2]);
//	System.out.println(VV[2][0] + " " + VV[2][1]+ " " + VV[2][2]);

	VV = MathTools.normalizeRows(VV);
//	System.out.println(VV[0][0] + " " + VV[0][1]+ " " + VV[0][2]);
//	System.out.println(VV[1][0] + " " + VV[1][1]+ " " + VV[1][2]);
//	System.out.println(VV[2][0] + " " + VV[2][1]+ " " + VV[2][2]);
		
	double normalizedCheck = MathTools.dot(VV[0],VV[0]);
	System.out.println("dot prod = " + normalizedCheck);

	double orthogCheck = MathTools.dot(VV[0],VV[11]);
	System.out.println("dot prod = " + orthogCheck);
	
	Matrix D = Eig.getD();
	double[] eigenvalue0 = D.transpose().getArray()[0];
	double[] eigenvalue1 = D.transpose().getArray()[1];
	double[] eigenvalue2 = D.transpose().getArray()[2];
	System.out.println(eigenvalue0[0] + " " + eigenvalue0[1]+ " " + eigenvalue0[2]);
	System.out.println(eigenvalue1[0] + " " + eigenvalue1[1]+ " " + eigenvalue1[2]);
	System.out.println(eigenvalue2[0] + " " + eigenvalue2[1]+ " " + eigenvalue2[2]);

	double[] e = Eig.getRealEigenvalues();
	System.out.println(e[0] + " " + e[1] + " " + e[2]);
    }
}

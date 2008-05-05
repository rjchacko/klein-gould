package rachele.ising.testCode.apps;

import scikit.numerics.Jama.*;
import rachele.util.MathTools;

public class jamaTestApp{

    public static void main(String[] argv){   
    int L = 3;
   	double [][] M = new double[L][L];
	
   	M[0][0] = 1;
   	M[0][1] = 2;
   	M[0][2] = 0;
   	M[1][0] = 0;
   	M[1][1] = 3;
   	M[1][2] = 0;
   	M[2][0] = 2;
   	M[2][1] = -4;
   	M[2][2] = 2;
//   	for(int x = 0; x < L; x++){
//   		for(int y = 0; y <= x; y++){
//   			M[x][y] = M[y][x] = x*y;
//   		}
//   	}

   	System.out.println("Here comes the matrix:");
   	for(int x = 0; x < L; x++){
   		for(int y = 0; y < L; y++){
   			System.out.print(" " + M[x][y] + " ");
   		}
   		System.out.println(" ");
   	}
	System.out.println(" ");
		
	Matrix A = new Matrix(M);	
	EigenvalueDecomposition Eig = A.eig();
	
	
	Matrix V = Eig.getV();
	double [][] VV = V.transpose().getArray();//take transpose to get eigenvectors as row

	System.out.println("Here comes V");
		for(int x = 0; x < L; x++){
     		for(int y = 0; y < L; y++){
     			System.out.print(" " + VV[x][y] + " ");
     		}
     		System.out.println(" ");
     	}
   		System.out.println(" ");
   		
	VV = MathTools.normalizeRows(VV);
	System.out.println("Here comes V normalized");
	for(int x = 0; x < L; x++){
 		for(int y = 0; y < L; y++){
 			System.out.print(" " + VV[x][y] + " ");
 		}
 		System.out.println(" ");
 	}
		System.out.println(" ");
		
	double normalizedCheck = MathTools.dot(VV[0],VV[0]);
	System.out.println("dot prod = " + normalizedCheck);

	double orthogCheck = MathTools.dot(VV[0],VV[1]);
	System.out.println("dot prod = " + orthogCheck);
	
	Matrix D = Eig.getD();
	double [][] DD = D.transpose().getArray();//take transpose to get eigenvectors as row
	System.out.println("Here comes D");
	for(int x = 0; x < L; x++){
 		for(int y = 0; y < L; y++){
 			System.out.print(" " + DD[x][y] + " ");
 		}
 		System.out.println(" ");
 	}
		System.out.println(" ");
	double[] eigenvalue0 = D.transpose().getArray()[0];
	double[] eigenvalue1 = D.transpose().getArray()[1];
	double[] eigenvalue2 = D.transpose().getArray()[2];
	System.out.println(eigenvalue0[0] + " " + eigenvalue0[1]+ " " + eigenvalue0[2]);
	System.out.println(eigenvalue1[0] + " " + eigenvalue1[1]+ " " + eigenvalue1[2]);
	System.out.println(eigenvalue2[0] + " " + eigenvalue2[1]+ " " + eigenvalue2[2]);

	double[] e = Eig.getRealEigenvalues();
	for(int i = 0; i <L; i++){
		System.out.println("ev " + i + " = " + e[i]);
	}
    
    //SVD
	SingularValueDecomposition Svd = A.svd();
	double [] sv = new double [L];
	sv = Svd.getSingularValues();

		for(int i = 0; i <L; i++){
			System.out.println("sv " + i + " = " + sv[i]);
   		}

		
		double sum = 0;
		int testInt=1;
	    for(int i = 0; i < L; i++){
	    		sum += M[testInt][i]*VV[testInt][i];
	    }     
	    double lambda = sum/VV[testInt][testInt];
	    System.out.println("ev = " + e[testInt] + " lambda = " + lambda);
    }

    
}

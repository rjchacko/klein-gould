package rachele.ising.dim2;

import rachele.util.FileUtil;
import rachele.util.FourierTransformer;
import rachele.util.MathTools;
import scikit.numerics.Jama.EigenvalueDecomposition;
import scikit.numerics.Jama.Matrix;
import scikit.numerics.fn.Function1D;
import scikit.numerics.fn.Function2D;
import scikit.util.DoubleArray;

public class StripeClumpFieldSim {
	int Lp;
	IsingField2D ising;
	FourierTransformer fft;
	int ky;
    double [] rhs, rhs2D, etaLT, eta, etaLT_k_slice, etaK; //right hand side
    double [] g_k, f_k, f_x, phi0_bar, c, eigenvalue;
	public double [] phi0;
    double [][] M;// Main matrix 
    public double [][] VV;// Eigenvalue matrix
	
	public StripeClumpFieldSim(IsingField2D ising, int ky){
		this.ky = ky;
		this.ising = ising;
		Lp = ising.Lp;
		fft = new FourierTransformer(Lp);

		g_k = new double[Lp];
		f_k = new double[Lp];
		f_x = new double [Lp];
		phi0 = new double [Lp];
		phi0_bar = new double [Lp];
		c = new double[Lp];
		eigenvalue = new double[Lp];
		M = new double [Lp][Lp];
		VV = new double [Lp][Lp];
		
		findPhi0andPhi0_bar();
		findMatrix();
		diagonalize();
		findAndDiagonalizeLittleMatrix();
	}

	public void findPhi0andPhi0_bar(){
		String fileName = "../../../research/javaData/configs1d/config";
		//need to make phi0 symmetric
		double [] tempPhi0 = FileUtil.readConfigFromFile(fileName, Lp);
		double minPhi0Value = 1.0;
		int minPhi0Location = -1;
		for (int i = 0; i < Lp; i++){
			if (tempPhi0[i] < minPhi0Value){
				minPhi0Location = i;
				minPhi0Value = tempPhi0[i];
				System.out.println(tempPhi0[i] + " " + i);
			}
		}	
		System.out.println(tempPhi0[minPhi0Location] + " " + minPhi0Location);
		for (int i = 0; i < Lp; i++){
			phi0[i] = tempPhi0[(minPhi0Location+i)%Lp];
			System.out.println("phi0 " + i + " = " + phi0[i]);
		}
		phi0_bar = fft.backConvolve1DwithFunction(phi0, new Function1D(){
			public double eval(double k1) {
				double kRx = 2*Math.PI*ising.R*k1/ising.L;
				double kRy = 0.0;
				if(ising.circleInteraction) return ising.findVkCircle(kRx*kRx + kRy*kRy);
				else return ising.findVkSquare(kRx, kRy);
			}
		});
		for(int i = 0; i < Lp; i++)
			f_x[i] = 1.0 / (1.0 - Math.pow(phi0[i],2));
		f_k = fft.calculate1DFT(f_x);
	}
	
	
	/**
	* Evaluates the linear theory prediction of 
	* configuration in Fourier space using a linear
	* combination of eigenvectors weighted by proper
	* coefficients and exp growth factors.
	*/
	public void etaLinearCombo(){
		for(int i = 0; i < Lp; i++){
			double sum = 0;
			for(int j = 0; j < Lp; j++){
				sum += c[j]*VV[j][i]*Math.exp(ising.time()*eigenvalue[j]);
			}
			//etaLC[i] = sum/(double)(Lp);
		}
	}
	
	/**
	* Evaluates the matrix, M, of the linear theory for the unstable 
	* (unmodified) Ising dynamics.
	*/
	public void findMatrix() {
		double kyValue = 2.0*Math.PI*ising.R*ky/ising.L;
		double kxValue;
		for (int i = 0; i < Lp; i++){
			for (int j = 0; j <= i; j++){
				M[i][j]=M[j][i]=-ising.T*f_k[(i-j+Lp)%Lp]/Lp;
			}
		}
		for (int i = 0; i < Lp; i++){
			if(i >= Lp/2)
				kxValue = 2.0*Math.PI*ising.R*(i-Lp)/ising.L;
			else
				kxValue = 2.0*Math.PI*ising.R*i/ising.L;
			if(ising.circleInteraction)
				M[i][i]-=ising.findVkCircle(Math.sqrt(kxValue*kxValue + kyValue*kyValue));
			else
				M[i][i]-=ising.findVkSquare(kxValue, kyValue);
		}
	}


	
	public void findAndDiagonalizeLittleMatrix() {
		
		int s = 4;//must be even
		double [][] MM = new double [s][s]; 
		
		for (int x = 0; x < s/2; x++){
			for (int y = 0; y < s/2; y++){
				MM[x][y]=M[x][y];
				MM[x+s/2][y]=M[x+Lp-s/2][y];
				MM[x][y+s/2]=M[x][y+Lp-s/2];
				MM[x+s/2][y+s/2]=M[x+Lp-s/2][y+Lp-s/2];
			}
		}

		Matrix mmatrix = new Matrix(MM);
		EigenvalueDecomposition Eigen;
		Eigen = mmatrix.eig();
		
		double [] eigenv = new double [s];
		eigenv = Eigen.getRealEigenvalues();
		

//		for (int i = 0; i < Lp; i ++)
//			System.out.println("eigenvalue " + i + " = " + eigenvalue[i]);
//		if(DoubleArray.max(eigenvalue)>0.0)
		System.out.println("eigenvalue little max ="  +  DoubleArray.max(eigenv));
		System.out.println("eigenvalue min ="  + DoubleArray.min(eigenv));

//		for(int i = Lp-1; i > Lp-5; i--){
//			for(int j = Lp-1; j > Lp-5;j--){
//				System.out.println(i + " " + j + " " + M[i][j] + " " + MM[i][j]);
//			}
//		}
	    
		
	}
	
//	public double [] simulateLinearK(){
//		double [] linearTheoryGrowth = new double[Lp*Lp];
//				
//		for (int i = 0; i < Lp; i++){
//			for (int j = 0; j < Lp; j++){
//				linearTheoryGrowth [i] += ising.dt*M[i][j]*etaLT_k_slice[j];
//			}
//		}
// 
//		return linearTheoryGrowth;
//	}

//	public double [] simulateLinearKbar(){
//		double kyValue = 2.0*Math.PI*ising.R*ky/ising.L;
//		double [] linearTheoryGrowth = new double[Lp*Lp];
//		double kxValue;
//		double [] f_bar = new double [Lp];
//		f_bar = fft.backConvolve1D(f_k, etaLT_k_slice);
//		convolve.registerLines("convolvution", new PointSet(1,1, f_bar), Color.green);
//		for (int i = 0; i < Lp; i++){
//			
//			if(i >= Lp/2)
//				kxValue = 2.0*Math.PI*ising.R*(i-Lp)/ising.L;
//			else
//				kxValue = 2.0*Math.PI*ising.R*i/ising.L;
//				
//			linearTheoryGrowth[i] = ising.dt*(-ising.findVkSquare(kxValue, kyValue)*etaLT_k_slice[i]);
//			linearTheoryGrowth [i] -= ising.dt*(ising.T*f_bar[i]);
//		}
//		etaLTkAcc.accum(ising.time(),Math.pow(etaLT_k_slice[0],2));		
//   		etaLTkAcc2.accum(ising.time(), Math.pow(etaLT_k_slice[1],2));
//   		etaLTkAcc3.accum(ising.time(),Math.pow(etaLT_k_slice[2],2));
//   		etaLTkAcc4.accum(ising.time(), Math.pow(etaLT_k_slice[3],2)); 
//		return linearTheoryGrowth;
//	}

	/**
	* Diagonalizes the matrix M.
	* Also some code that can be used to check for
	* orthonormality of eigenvectors. (Eigenvectors should
	* be orthonormal for Hermitian matrices.  Matrix for unmodified
	* dynamics is Hermitian, but matrix for modified dynamcis
	* is not.
	*/
	public void diagonalize(){
		Matrix matrix = new Matrix(M);
		EigenvalueDecomposition Eig;
		Eig = matrix.eig();
		eigenvalue = new double [Lp];
		eigenvalue = Eig.getRealEigenvalues();
		

		for (int i = 0; i < Lp; i ++)
			System.out.println("eigenvalue " + i + " = " + eigenvalue[i]);
		if(DoubleArray.max(eigenvalue)>0.0)System.out.println("eigenvalue max ="  +  DoubleArray.max(eigenvalue));
		System.out.println("eigenvalue min ="  + DoubleArray.min(eigenvalue));
		Matrix V = Eig.getV();
		VV = V.transpose().getArray();
		VV = MathTools.normalizeRows(VV);
		for (int i = 0; i < Lp; i++){
			double normalizedCheck = MathTools.dot(VV[i],VV[i]);
			if(Math.abs(normalizedCheck-1.0)>.000001)
				System.out.println("dot prod norm = " + normalizedCheck);
			for (int j = 0; j < Lp; j++){
				double orthogCheck = MathTools.dot(VV[i],VV[j]);
				if(i != j & Math.abs(orthogCheck)>0.000001)
					System.out.println("dot prod orth = " + orthogCheck);			
			}
		}

		//find coeffieicnts
		//normalize the current eta slice
//		double [] normedEta = MathTools.normalize(etaLT_k_slice);
//		double normCheck = 0;
//		for (int i = 0; i < Lp; i ++){
//			c[i] = MathTools.dot(normedEta, VV[i]);
//			normCheck += c[i]*c[i];
//		}
//		
//		double sum = 0;
//		int testInt=1;
//	    for(int i = 0; i < Lp; i++){
//	    		sum += M[testInt][i]*VV[testInt][i];
//	    } 
//	    double lambda = sum/VV[testInt][testInt];
//	    System.out.println("ev = " + eigenvalue[testInt] + " lambda = " + lambda);
	    
	}
	public double [] simulateLinearMod_kSpace(int ky){
		double [] change = new double [Lp*Lp];

//		for (int ky = -Lp/2; ky < Lp/2; ky++){
//		double kyValue = 2.0*Math.PI*ising.R*ky/ising.L;
			for(int kx = -Lp/2; kx < Lp/2; kx++){	
				//int i = Lp*((ky+Lp)%Lp) + (kx+Lp)%Lp;
				double sum = 0;				
				for (int kxx = -Lp/2; kxx < Lp/2; kxx++){
					//double kxxValue = 2.0*Math.PI*ising.R*kxx/ising.L;
					//int eta_k_int = Lp*((ky+Lp)%Lp) + (kxx+Lp)%Lp;
					sum += M[(kx+Lp)%Lp][(kxx+Lp)%Lp]*etaLT_k_slice[(kxx+Lp)%Lp];				
				}
				change[(kx+Lp)%Lp] = sum/Lp;
			}
	//	}
		return change;
	}
	
	
	public double [] simulateLinear(double [] etaLinear){
		double [] linearTheoryGrowth = new double[Lp*Lp];
		double [] eta_bar = new double[Lp*Lp];
		double [] phi0_2D = new double[Lp*Lp];
		for(int i = 0; i < Lp*Lp; i++)
			phi0_2D[i] = phi0[i%Lp];
//		ising.convolveWithRange(eta, eta_bar, ising.R);
		eta_bar = fft.convolve2DwithFunction(etaLinear, new Function2D(){
			public double eval(double k1, double k2) {
				double kRx = 2*Math.PI*ising.R*k1/ising.L;
				double kRy = 2*Math.PI*ising.R*k2/ising.L;
				if(ising.circleInteraction) return ising.findVkCircle(kRx*kRx + kRy*kRy);
				else return ising.findVkSquare(kRx, kRy);
			}
		});

		for (int i = 0; i < Lp*Lp; i++)
			linearTheoryGrowth[i] = ising.dt*(-eta_bar[i]-ising.T*etaLinear[i]*f_x[i%Lp]);
		return linearTheoryGrowth;
	}
	
	
//	/**
//	* Evaluates the matrix, M, of the linear theory for the stable 
//	* (modified) Ising dynamics.  Not working at present.
//	*/
//	public void findModMatrix(){
////		double kyValue = 2.0*Math.PI*ising.R*ky/ising.L;
////		double kxValue;
////		for (int i = 0; i < Lp; i++){
////			for (int j = 0; j <= i; j++){
////				M[i][j]=M[j][i]=-ising.T*f_k[(i-j+Lp)%Lp]/Lp;
////			}
////		}
////		for (int i = 0; i < Lp; i++){
////			if(i >= Lp/2)
////				kxValue = 2.0*Math.PI*ising.R*(i-Lp)/ising.L;
////			else
////				kxValue = 2.0*Math.PI*ising.R*i/ising.L;
////			M[i][i]-=ising.findVkSquare(kxValue, kyValue);
////		}
////		
//		for (int i = 0; i < Lp; i++){
//			double alpha_x = 1.0-Math.pow(phi0[i],2);
//			G[i] = -alpha_x*alpha_x;
//			F[i] = alpha_x*(-ising.T-4*phi0[i]*(-phi0_bar[i]-ising.T*scikit.numerics.Math2.atanh(phi0[i]) + ising.H));				
//		}
////		g_k = fft.calculate1DFT(g);
////		f_k = fft.calculate1DFT(f);
////
////		//for (int ky = -Lp/2; ky < Lp/2; ky++){
////		//double kyValue = 2.0*Math.PI*ising.R*ky/ising.L;
////		for(int kx = -Lp/2; kx < Lp/2; kx++){	
////			//int i = Lp*((ky+Lp)%Lp) + (kx+Lp)%Lp;
////			//double sum = 0;				
////			for (int kxx = -Lp/2; kxx < Lp/2; kxx++){
////				double kxxValue = 2.0*Math.PI*ising.R*kxx/ising.L;
////				//int eta_k_int = Lp*((ky+Lp)%Lp) + (kxx+Lp)%Lp;
////				M[(kx+Lp)%Lp][(kxx+Lp)%Lp] += (g_k [(kxx-kx+Lp)%Lp] * ising.findVkSquare(kxxValue,kyValue) + f_k[(kxx-kx+Lp)%Lp]);				
////			}
////		}
//	}	
}

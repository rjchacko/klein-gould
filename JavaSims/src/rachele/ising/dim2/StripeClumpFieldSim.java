package rachele.ising.dim2;

//import java.io.File;

//import java.awt.Color;
import java.io.File;
//import java.security.Policy.Parameters;

import rachele.util.FileUtil;
import rachele.util.FourierTransformer;
import rachele.util.MathTools;
//import scikit.dataset.PointSet;
import scikit.jobs.params.Parameters;
import scikit.numerics.Jama.EigenvalueDecomposition;
import scikit.numerics.Jama.Matrix;
import scikit.numerics.fn.Function1D;
import scikit.numerics.fn.Function2D;
import scikit.util.DoubleArray; 

/**
* Calculates matrix, eigenvectors, eigenvalues....
* 
* The simulation part of this simulation has not
* 
*/
public class StripeClumpFieldSim {
	int Lp;
	IsingField2D ising;
	FourierTransformer fft;
	int ky;
    double [] rhs, rhs2D,  etaLT_k_slice, etaK; //right hand side
    double [] g_k, f_x, phi0_bar, c, f_G;
	public double [] phi0, f_k, f_Gk, eigenvalue, etaLT, etaLT2, etaLT_k, etaLT2D_k, etaKchange, etaKchange2;
    public double [] etaBar_k, etaBarCheck, eta_bar2;
	double [][] M;// Main matrix 
    public double [][] VV;// Eigenvalue matrix
	String phi0file;
    
	public StripeClumpFieldSim(IsingField2D ising, Parameters params){
		phi0file = params.sget("1D Input File");
		ky = params.iget("ky");
		this.ising = ising;
		Lp = ising.Lp;
		fft = new FourierTransformer(Lp);

		g_k = new double[Lp];
		f_Gk = new double[Lp];
		f_G = new double[Lp];
		f_k = new double[Lp];
		f_x = new double [Lp];
		etaKchange = new double [Lp*Lp];
		etaKchange2 = new double [Lp*Lp];
		phi0 = new double [Lp];
		phi0_bar = new double [Lp];
		c = new double[Lp];
		eigenvalue = new double[Lp];
		M = new double [Lp][Lp];
		VV = new double [Lp][Lp];
		
		findPhi0andPhi0_bar();
		if (params.sget("Dynamics")=="Glauber"){
			calcXtraGlauberFunctions();
			findGlauberMatrix();
		}
		else if (params.sget("Dynamics")=="Langevin"){
			calcXtraFunctions();			
			findMatrix();
		}
		diagonalize();
		for (int i=0; i < Lp; i ++){
			System.out.println("eigenvalue " + i + " = " + eigenvalue[i]);
		}
		
		
		writeMatrixResultsToFile(params);

		//findAndDiagonalizeLittleMatrix();
	}

	public void writeMatrixResultsToFile(Parameters params){
		String outDir = params.sget("Data Dir");
		StringBuffer sb = new StringBuffer();
		sb.append(outDir);	sb.append(File.separator); sb.append("L"); sb.append(ising.Lp);
		int range = (int)(ising.Lp/params.fget("L/R")); sb.append("R"); sb.append(range);
		sb.append("T"); sb.append(ising.T); sb.append("h"); sb.append(ising.H);
		String outFilePreamble = sb.toString();
		String evalueFile = outFilePreamble + "eValues";	
		StringBuffer mb = new StringBuffer();
		mb.append("# kR line value is ");
		double kRLine = 2*ising.R*Math.PI*(params.iget("ky"))/ising.L; mb.append(kRLine);
		FileUtil.initFile(evalueFile, params, mb.toString());
		FileUtil.printlnToFile(evalueFile, "# Results of StripeClumpFieldSim calculations ");
		FileUtil.printlnToFile(evalueFile, "inputFile = " + phi0file);
		for (int i=0; i < Lp; i ++){
			//System.out.println("start " + i);
			FileUtil.printlnToFile(evalueFile, i, eigenvalue[i]);
		}

		double ratio = eigenvalue[Lp-1]/eigenvalue[Lp-2];
		FileUtil.printlnToFile(evalueFile, "# Ratio of highest 2 eigenvalues = " + ratio);
		int eVectorNo = 5;
		for (int i = 0; i < eVectorNo; i++){
			String evectorFile = outFilePreamble + "eVec" + i;
			FileUtil.initFile(evectorFile, params);
			int no = Lp-i-1;
			FileUtil.printlnToFile(evectorFile, "# eigenvector no ", no);
			for (int j = 0; j < Lp; j++)
				FileUtil.printlnToFile(evectorFile, j, VV[no][j]);				
			sb.deleteCharAt(sb.length()-1);
		}
		for (int i=0; i < Lp; i ++){
			System.out.println("eigenvalue " + i + " = " + eigenvalue[i]);
		}
		System.out.println("kR line = " + kRLine);


	}
	
	public void findPhi0andPhi0_bar(){
		//String fileName = "../../../research/javaData/configs1d/config";
		//need to make phi0 symmetric
		double [] tempPhi0 = FileUtil.readConfigFromFile(phi0file, Lp);
		double minPhi0Value = 1.0;
		int minPhi0Location = -1;
		for (int i = 0; i < Lp; i++){
			if (tempPhi0[i] < minPhi0Value){
				minPhi0Location = i;
				minPhi0Value = tempPhi0[i];
				//System.out.println(tempPhi0[i] + " " + i);
			}
		}	
		//System.out.println(tempPhi0[minPhi0Location] + " " + minPhi0Location);
		for (int i = 0; i < Lp; i++){
			phi0[i] = tempPhi0[(minPhi0Location+i)%Lp];
			//System.out.println("phi0 " + i + " = " + phi0[i]);
		}
		phi0_bar = fft.backConvolve1DwithFunction(phi0, new Function1D(){
			public double eval(double k1) {
				double kRx = 2*Math.PI*ising.R*k1/ising.L;
				double kRy = 0.0;
				if(ising.circleInteraction) return ising.findVkCircle(kRx*kRx + kRy*kRy);
				else return ising.findVkSquare(kRx, kRy);
			}
		});
		

	}
	
	public void calcXtraFunctions(){
		for(int i = 0; i < Lp; i++)
			f_x[i] = 1.0 / (1.0 - Math.pow(phi0[i],2));
		f_k = fft.calculate1DFT(f_x);		
	}
	
	public void calcXtraGlauberFunctions(){
		for(int i = 0; i < Lp; i++){
//			double tanhTerm = Math.tanh(phi0_bar[i]/ising.T+ising.H/ising.T);
//			f_G[i] = tanhTerm*tanhTerm;
			f_G[i] = phi0[i]*phi0[i];
//			System.out.println("tt " +tanhTerm);
		}
		f_Gk = fft.calculate1DFT(f_G);
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
	public void findMatrix_mobility1() {
		double kyValue = 2.0*Math.PI*ising.R*ky/ising.L;
		double kxValue;
		for (int i = 0; i < Lp; i++){
			for (int j = 0; j <= i; j++){
				M[i][j] = M[j][i] = (-ising.T*f_k[(i-j+Lp)%Lp]);
			}
		}
		for (int i = 0; i < Lp; i++){
			if(i >= Lp/2)
				kxValue = 2.0*Math.PI*ising.R*(i-Lp)/ising.L;
			else
				kxValue = 2.0*Math.PI*ising.R*i/ising.L;
			if(ising.circleInteraction)
				M[i][i] += ising.J*ising.findVkCircle(Math.sqrt(kxValue*kxValue + kyValue*kyValue));
			else
				M[i][i] += ising.J*ising.findVkSquare(kxValue, kyValue);
		}
	}
	
//	public void findMatrix() {
//		double kyValue = 2.0*Math.PI*ising.R*ky/ising.L;
//		double kxValue;
//		double [] mobility1d = new double [Lp];
//		double [] f2 = new double [Lp];
//		for (int i = 0; i < Lp; i++){
//			mobility1d[i] = ising.mobility[i+Lp*ky];
//			f2[i] = mobility1d[i]*ising.T/(1-phi0[i]*phi0[i]);
//		}
//		double [] mobility_k = fft.calculate1DFT(mobility1d);
//		double [] f2k = fft.calculate1DFT(f2);
//
//		for (int i = 0; i < Lp; i++){
//			if(i >= Lp/2)
//				kxValue = 2.0*Math.PI*ising.R*(i-Lp)/ising.L;
//			else
//				kxValue = 2.0*Math.PI*ising.R*i/ising.L;
//			for (int j = 0; j <= i; j++){
////				M[i][j] = M[j][i] = mobility_k[(i-j+Lp)%Lp]*ising.J*ising.findVkSquare(kxValue, kyValue);
////				M[i][j] = M[j][i] = mobility_k[(j-i+Lp)%Lp]*ising.J*ising.findVkSquare(kxValue, kyValue)- f2k[(j-i+Lp)%Lp];
//				M[i][j] = M[j][i] = - f2k[(j-i+Lp)%Lp];
//
//			}
////			M[i][i] += mobility_k[0]*ising.J*ising.findVkSquare(kxValue, kyValue);
//		}
//
//	}
	
	public void findMatrix() {
		//This is final and working.  Don't fuck it up.
		double [] mobility1d = new double [Lp];
		double [] f2 = new double [Lp];
		for (int i = 0; i < Lp; i++){
			mobility1d[i] = ising.mobility[i+Lp*ky];
			f2[i] = mobility1d[i]*ising.T/(1-phi0[i]*phi0[i]);
		}
		double [] mobility_k = fft.calculate1DFT(mobility1d);
		double [] f2k = fft.calculate1DFT(f2);

		int y1=ky;
		for (int x1 = -Lp/2; x1 < Lp/2; x1++) {
			for (int x2 = -Lp/2; x2 < Lp/2; x2++) {
				double kRx = 2*Math.PI*ising.R*x2/ising.L;
				double kRy = 2*Math.PI*ising.R*y1/ising.L;
				M[(x2+Lp)%Lp][(x1+Lp)%Lp] = mobility_k[(x1-x2+Lp)%Lp]*ising.J*ising.findVkSquare(kRx, kRy) - f2k[(x1-x2+Lp)%Lp];
			}
		}

	}

	
	public void findGlauberMatrix() {
		double kyValue = 2.0*Math.PI*ising.R*ky/ising.L;
		double kxValue;
		for (int i = 0; i < Lp; i++){
			if(i >= Lp/2)
				kxValue = 2.0*Math.PI*ising.R*(i-Lp)/ising.L;
			else
				kxValue = 2.0*Math.PI*ising.R*i/ising.L;
			for (int j = 0; j <= i; j++)
				M[i][j] = M[j][i] = - ising.J*ising.findVkSquare(kxValue, kyValue)*f_Gk[(j-i+Lp)%Lp]/(ising.T);
		}
		for (int i = 0; i < Lp; i++){
			if(i >= Lp/2)
				kxValue = 2.0*Math.PI*ising.R*(i-Lp)/ising.L;
			else
				kxValue = 2.0*Math.PI*ising.R*i/ising.L;
			M[i][i] += -1 + ising.J*ising.findVkSquare(kxValue, kyValue)/ising.T;
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
		System.out.println("eigenvalue little max ="  +  DoubleArray.max(eigenv));
		System.out.println("eigenvalue min ="  + DoubleArray.min(eigenv));
	}
	
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
			//System.out.println("eigenvalue " + i + " = " + eigenvalue[i]);
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
	    
	}
	
	public void initEta(){
		etaLT = new double [Lp*Lp];
		etaLT2 = new double [Lp*Lp];
		for (int i = 0; i < Lp*Lp; i++)
			etaLT[i] = etaLT2[i] = ising.phi[i] - phi0[i%Lp];
		etaLT2D_k = new double [Lp*Lp];
		etaLT2D_k = fft.calculate2DFT(etaLT);
		etaLT_k = new double [Lp];
		etaBar_k = new double [Lp*Lp];
		etaBarCheck = new double [Lp*Lp];
		eta_bar2 = new double[Lp*Lp];

		
//		for (int i = ky*Lp; i < Lp*(ky+1); i++)
//			etaLT_k[i-ky*Lp] = etaLT2D_k[i]; 
		for (int i = 0; i < Lp; i++)
			etaLT_k[i] = etaLT2D_k[ky*Lp+i];
		System.out.println("init eta");
	}


	/**
	 * 
	 */
	public void simulateLinear(){
		double [] linearTheoryGrowth = new double[Lp*Lp];
		double [] eta_bar = new double[Lp*Lp];
		double [] phi0_2D = new double[Lp*Lp];
		for(int i = 0; i < Lp*Lp; i++)
			phi0_2D[i] = phi0[i%Lp];
		eta_bar = fft.convolve2DwithFunction(etaLT, new Function2D(){
			public double eval(double k1, double k2) {
				double kRx = 2*Math.PI*ising.R*k1/ising.L;
				double kRy = 2*Math.PI*ising.R*k2/ising.L;
				if(ising.circleInteraction) return ising.findVkCircle(kRx*kRx + kRy*kRy);
				else return ising.findVkSquare(kRx, kRy);
			}
		});
		for (int i = 0; i < Lp*Lp; i++){
//			linearTheoryGrowth[i] = ising.dt*(-etaLT[i]);
			linearTheoryGrowth[i] = ising.dt*ising.mobility[i]*(-ising.T*etaLT[i]*f_x[i%Lp]);
			linearTheoryGrowth[i] += ising.dt*(-eta_bar[i]*ising.mobility[i]);
		}
		

		
		for (int i = 0; i < Lp*Lp; i++)
			etaLT[i] += linearTheoryGrowth[i];

//		//this works:
//		double [] linearGk = fft.calculate2DFT(linearTheoryGrowth);
//		for (int i = 0; i < Lp*Lp; i++){
//				etaLT2D_k[i] += linearGk[i];	
//		}
//
//		//this works
//		for (int i = 0; i < Lp; i++){
//			etaLT_k[i] += linearGk[i+Lp*ky];	
//		}
//		
		

			
	}
	
	public void simulateLinearK(){
		double [] linearTheoryGrowth = new double[Lp];
		for (int i = 0; i < Lp; i++){
			for (int j = 0; j < Lp; j++)
				linearTheoryGrowth[i] += ising.dt*(M[j][i]*etaLT_k[j]);
		}
		for (int i = 0; i < Lp; i++){
			etaLT_k[i] += linearTheoryGrowth[i];
		}	
	}
	
	public void simulateLinearKbar1(){
		//this works
		
		double [] eta_bar = new double[Lp*Lp];
		double [] phi0_2D = new double[Lp*Lp];
		for(int i = 0; i < Lp*Lp; i++)
			phi0_2D[i] = phi0[i%Lp];
		double [] mobility1d = new double [Lp];
		for(int i = 0; i <Lp; i++)
			mobility1d[i] = ising.mobility[i];
		double [] mobility1dk = fft.calculate1DFT(mobility1d);
		eta_bar = fft.convolve2DwithFunction(etaLT, new Function2D(){
			public double eval(double k1, double k2) {
				double kRx = 2*Math.PI*ising.R*k1/ising.L;
				double kRy = 2*Math.PI*ising.R*k2/ising.L;
				if(ising.circleInteraction) return ising.findVkCircle(kRx*kRx + kRy*kRy);
				else return ising.findVkSquare(kRx, kRy);
			}
		});
		double [] etaBar_k = fft.calculate2DFT(eta_bar);
		double [] etakDot = new double [Lp*Lp];
		

		for (int i = 0; i < Lp; i++){
			for (int j = 0; j < Lp; j++){
				etakDot[i] += -mobility1dk[(i-j+Lp)%Lp]*etaBar_k[j+Lp*ky];
			}
			etakDot[i] -=etaLT_k[i];
		}
		for (int i = 0; i < Lp; i++){
			etaLT_k[i] += ising.dt*etakDot[i];
			System.out.println(i + " " + etakDot[i]);
		}
	}

	public void simulateLinearKbar2d(){
		//this works
		
		double [] eta_bar = new double[Lp*Lp];
		double [] phi0_2D = new double[Lp*Lp];
		for(int i = 0; i < Lp*Lp; i++)
			phi0_2D[i] = phi0[i%Lp];
		double [] mobility1d = new double [Lp];
		for(int i = 0; i <Lp; i++)
			mobility1d[i] = ising.mobility[i];
		double [] mobility1dk = fft.calculate1DFT(mobility1d);
		eta_bar = fft.convolve2DwithFunction(etaLT, new Function2D(){
			public double eval(double k1, double k2) {
				double kRx = 2*Math.PI*ising.R*k1/ising.L;
				double kRy = 2*Math.PI*ising.R*k2/ising.L;
				if(ising.circleInteraction) return ising.findVkCircle(kRx*kRx + kRy*kRy);
				else return ising.findVkSquare(kRx, kRy);
			}
		});
		double [] etaBar_k = fft.calculate2DFT(eta_bar);
		double [] etakDot = new double [Lp*Lp];
		

		//but etabar k = vk*etak
		
		for (int i = 0; i < Lp*Lp; i++){
			int thisky = i/Lp;
			int thiskx = i%Lp;
			for (int j = 0; j < Lp; j++){
				etakDot[i] += -mobility1dk[(thiskx-j+Lp)%Lp]*etaBar_k[j+Lp*thisky];
			}
			etakDot[i] -=etaLT2D_k[i];
		}
		for (int i = 0; i < Lp*Lp; i++){
			etaLT2D_k[i] += ising.dt*etakDot[i];
			System.out.println(i + " " + etakDot[i]);
		}
	}


	public void simulateLinearKbar2d2(){
		//this works
		
//		double [] eta_bar = new double[Lp*Lp];
		double [] phi0_2D = new double[Lp*Lp];
		for(int i = 0; i < Lp*Lp; i++)
			phi0_2D[i] = phi0[i%Lp];
		double [] mobility1d = new double [Lp];
		for(int i = 0; i <Lp; i++)
			mobility1d[i] = ising.mobility[i];
		double [] mobility1dk = fft.calculate1DFT(mobility1d);
		eta_bar2 = fft.convolve2DwithFunction2(etaLT, new Function2D(){
			public double eval(double k1, double k2) {
				double kRx = 2*Math.PI*ising.R*k1/ising.L;
				double kRy = 2*Math.PI*ising.R*k2/ising.L;
				if(ising.circleInteraction) return ising.findVkCircle(kRx*kRx + kRy*kRy);
				else return ising.findVkSquare(kRx, kRy);
			}
		});
//		etaBar_k = new double[Lp*Lp];
		etaBar_k = fft.calculate2DFT(eta_bar2);
//		etaBarCheck = new double [Lp*Lp];
		double [] eta2dk= fft.calculate2DFT(etaLT);
		for (int y = -Lp/2; y < Lp/2; y++) {
			for (int x = -Lp/2; x < Lp/2; x++) {
				int i = (x+Lp)%Lp + Lp*((y+Lp)%Lp);
				double kRx = 2*Math.PI*ising.R*x/ising.L;
				double kRy = 2*Math.PI*ising.R*y/ising.L;
				etaBarCheck[i] = eta2dk[i]*ising.findVkSquare(kRx, kRy);				
			}
		}
//		etaBarCheck = fft.calculate2DBackFT(etaBarCheck);
		double [] etakDot = new double [Lp*Lp];
		


		//this works

		for (int y1 = -Lp/2; y1 < Lp/2; y1++) {
			for (int x1 = -Lp/2; x1 < Lp/2; x1++) {
				int i1 = (x1+Lp)%Lp + Lp*((y1+Lp)%Lp);
					for (int x2 = -Lp/2; x2 < Lp/2; x2++) {
						int i2 = (x2+Lp)%Lp + Lp*((y1+Lp)%Lp);
						double kRx = 2*Math.PI*ising.R*x2/ising.L;
						double kRy = 2*Math.PI*ising.R*y1/ising.L;
						etakDot[i1] += -mobility1dk[(x1-x2+Lp)%Lp]*eta2dk[i2]*ising.findVkSquare(kRx, kRy);
					}

				etakDot[i1] -=etaLT2D_k[i1];
			}
		}
		
	
		
//		for (int i = 0; i < Lp*Lp; i++){
//			int thisky = i/Lp;
//			int thiskx = i%Lp;
//			for (int j = 0; j < Lp; j++){
//				etakDot[i] += -mobility1dk[(thiskx-j+Lp)%Lp]*etaBarCheck[j+Lp*thisky];
//			}
//			etakDot[i] -=etaLT2D_k[i];
//		}
		for (int i = 0; i < Lp*Lp; i++){
			etaLT2D_k[i] += ising.dt*etakDot[i];

		}
	}	

	public void simulateLinearKbar1d(){
		//this works

		double [] phi0_2D = new double[Lp*Lp];
		for(int i = 0; i < Lp*Lp; i++)
			phi0_2D[i] = phi0[i%Lp];
		double [] mobility1d = new double [Lp];
		for(int i = 0; i <Lp; i++)
			mobility1d[i] = ising.mobility[i];
		double [] mobilityk = fft.calculate1DFT(mobility1d);
		eta_bar2 = fft.convolve2DwithFunction2(etaLT, new Function2D(){
			public double eval(double k1, double k2) {
				double kRx = 2*Math.PI*ising.R*k1/ising.L;
				double kRy = 2*Math.PI*ising.R*k2/ising.L;
				if(ising.circleInteraction) return ising.findVkCircle(kRx*kRx + kRy*kRy);
				else return ising.findVkSquare(kRx, kRy);
			}
		});
		etaBar_k = fft.calculate2DFT(eta_bar2);
		double [] eta2dk= fft.calculate2DFT(etaLT);
		for (int y = -Lp/2; y < Lp/2; y++) {
			for (int x = -Lp/2; x < Lp/2; x++) {
				int i = (x+Lp)%Lp + Lp*((y+Lp)%Lp);
				double kRx = 2*Math.PI*ising.R*x/ising.L;
				double kRy = 2*Math.PI*ising.R*y/ising.L;
				etaBarCheck[i] = eta2dk[i]*ising.findVkSquare(kRx, kRy);				
			}
		}
		double [] etakDot = new double [Lp];


		
		//this works
//		double [] f2 = new double [Lp];
//		for (int i = 0; i < Lp; i++){
//			mobility1d[i] = ising.mobility[i+Lp*ky];
//			f2[i] = mobility1d[i]*ising.T/(1-phi0[i]*phi0[i]);
//		}
//		double [] mobility_k = fft.calculate1DFT(mobility1d);
//		double [] f2k = fft.calculate1DFT(f2);
//				for (int i = 0; i < Lp; i++){
//			for (int j = 0; j < Lp; j++){
//			M[i][j] = M[j][i] = - f2k[(j-i+Lp)%Lp];
//			}
//		}
		
		
//		for (int x1 = -Lp/2; x1 < Lp/2; x1++) {
//			for (int x2 = -Lp/2; x2 < Lp/2; x2++) {
//				double kRx = 2*Math.PI*ising.R*x2/ising.L;
//				double kRy = 2*Math.PI*ising.R*ky/ising.L;
//				M[(x1+Lp)%Lp][(x2+Lp)%Lp] = M[(x2+Lp)%Lp][(x1+Lp)%Lp] = mobility_k[(x1-x2+Lp)%Lp]*ising.J*ising.findVkSquare(kRx, kRy);//-f2k[(x1-x2+Lp)%Lp];//



		
		int y1= ky;
		for (int x1 = -Lp/2; x1 < Lp/2; x1++) {
			for (int x2 = -Lp/2; x2 < Lp/2; x2++) {
				double kRx = 2*Math.PI*ising.R*x2/ising.L;
				double kRy = 2*Math.PI*ising.R*y1/ising.L;
//				etakDot[(x1+Lp)%Lp] += mobilityk[(x1-x2+Lp)%Lp]*ising.J*ising.findVkSquare(kRx, kRy)*etaLT_k[(x2+Lp)%Lp];
				M[(x1+Lp)%Lp][(x2+Lp)%Lp] = mobilityk[(x1-x2+Lp)%Lp]*ising.J*ising.findVkSquare(kRx, kRy);
			}


		}

		for (int x1 = -Lp/2; x1 < Lp/2; x1++) {
			M[(x1+Lp)%Lp][(x1+Lp)%Lp] -= 1.0;//etaLT_k[(x1+Lp)%Lp];
//			etakDot[(x1+Lp)%Lp] -=etaLT_k[(x1+Lp)%Lp];			
		}
		
		//This does work
//		for (int x1 = -Lp/2; x1 < Lp/2; x1++) {
//			for (int x2 = -Lp/2; x2 < Lp/2; x2++) {
//				etakDot[(x1+Lp)%Lp] += M[(x1+Lp)%Lp][(x2+Lp)%Lp]*etaLT_k[(x2+Lp)%Lp]; 
//			}
//		}
		
//		This is working now
		for (int i = 0; i < Lp; i++){
			for (int j = 0; j < Lp; j++){
				etakDot[i] += M[i][j]*etaLT_k[j];
			}
		}


		
		for (int i = 0; i < Lp; i++){
			etaLT_k[i] += ising.dt*etakDot[i];

		}
	}	

}

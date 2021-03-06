package rachele.ising.dim2.apps;

import static scikit.util.Utilities.asList;

import java.io.DataInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.awt.Color;
import rachele.util.FourierTransformer;
import rachele.util.MathTools;
import rachele.ising.dim2.IsingField2Dopt;
import rachele.util.FileUtil;
import scikit.graphics.dim2.Geom2D;
import scikit.graphics.dim2.Grid;
import scikit.graphics.dim2.Plot;
import scikit.jobs.Control;
import scikit.jobs.Job;
import scikit.jobs.Simulation;
import scikit.jobs.params.ChoiceValue;
import scikit.jobs.params.DoubleValue;
import scikit.numerics.Jama.*;
import scikit.numerics.fn.Function1D;
import scikit.numerics.fn.Function2D;
import scikit.util.DoubleArray;
import scikit.dataset.Accumulator;
import scikit.dataset.PointSet;
import kip.util.Random;

/**
* Compares the ising simulations to the linear theory
* following a quench in external field.
* 
* This works for the circle shaped interaction shape.  
* Also can use TestLinearApp.java with fixed length scales 
* for square shape interaction.
* 
* The linear theory is still:
* \eta_dot(x,y) = -\int dx' V(x-x',y-y') \eta(x',y')
* 				-T \eta(x,y)/(1-phi_0(x)^2)
* Unlike with the square interaction shape,
* the circle interaction does not break down neatly into slices 
* in Fourier space.
* It is still a matrix equation. in Fourier space:  
* \etaVec[i] = M[i,j]\etaVec[j]
* where \etaVec is a vector of length Lp*Lp; each dimension
* represents an allowed value of (kx, ky).
*/
public class TestLinearOptApp extends Simulation{
    IsingField2Dopt ising;
    FourierTransformer fft;
    double [] rhs, rhs2D, etaLT, eta, etaLT_k_slice; //right hand side
    double [] g_k, c, eigenvalue, etaLC;
    double [] phi0, phi0_bar; // Background stripe configuration and this configuration convoluted with potential.
    double [] f_x, f_k; // f(x) = 1/(1-phi_0^2)
    double [] F, G; //Matrices for modified dynamics
    double [][] M, VV; //Main matrix and eigenvalue matrix
    int ky;
    double kRChunk; //=2piR/L
    boolean clearFile;
	Random random = new Random();
    
    //RUN OPTIONS
    boolean accEtaValues = false;
    boolean modifiedDynamics = false;
    boolean writeToFile = true;
    
    public Accumulator etaAcc;
    public Accumulator etaAcc2;
    public Accumulator etaAcc3;
    public Accumulator etaAcc4;

    public Accumulator etaLTAcc; 
    public Accumulator etaLTAcc2;
    public Accumulator etaLTAcc3;
    public Accumulator etaLTAcc4;
    
    public Accumulator etaLTkAcc; 
    public Accumulator etaLTkAcc2;
    public Accumulator etaLTkAcc3;
    public Accumulator etaLTkAcc4;

    public Accumulator etaLCAcc; 
    public Accumulator etaLCAcc2;
    public Accumulator etaLCAcc3;
    public Accumulator etaLCAcc4;
    
    public int Lp;
    Grid etaDot = new Grid("ising delta phi");
    Grid eta_k_LTGrid = new Grid("rhs");
    Grid phiGrid = new Grid("phi field");
    Grid etaDotSF = new Grid("sf phi field");
    Plot etakPlot = new Plot("eta(k)");
    Grid phi0Grid = new Grid("phi_0");
    Grid phi0SliceGrid = new Grid("phi_0 slice");
	Plot hSlice = new Plot("Horizontal Slice");
	Plot vSlice = new Plot("Vertical Slice"); 
    Plot etaVsTimeSim = new Plot("eta v t");
    Plot etaVsTimeLinearK = new Plot("eta v t LT k space");
    Plot etaVsTimeLinear = new Plot("eta v t LT");
    Plot etaVsTimeLC = new Plot("eta v t LC");
    Plot fkPlot = new Plot ("f(k)");
    Plot etakSimPlot = new Plot("eta k sim");
    Plot convolve = new Plot("convolve");

    
	public static void main(String[] args) {
		new Control(new TestLinearOptApp(), "Ising Linear Test");
	}
	public void load(Control c) {

		if(accEtaValues) c.frameTogether("accs", etaVsTimeSim, etaVsTimeLinear, etaVsTimeLinearK, etaVsTimeLC);
		c.frame(phiGrid);
		c.frameTogether("Slices", vSlice, hSlice);
		params.addm("Zoom", new ChoiceValue("Yes", "No"));
		params.addm("Interaction", new ChoiceValue( "Circle","Square" ));
		params.addm("Theory", new ChoiceValue("Exact", "Slow Near Edge", "Dynamic dt"));
		params.addm("Dynamics?", new ChoiceValue("Langevin No M Convervation", "Langevin Conserve M"));
		params.add("Init Conditions", new ChoiceValue("Read in 1D solution","Read From File", "Random Gaussian" ));
		params.addm("Noise", new DoubleValue(0, 0, 1.0).withSlider());
		params.addm("Horizontal Slice", new DoubleValue(0.5, 0, 0.9999).withSlider());
		params.addm("Vertical Slice", new DoubleValue(0.5, 0, 0.9999).withSlider());
		params.addm("T", 0.06);
		params.addm("H", 0.5);
		params.addm("Rx", 2490000.0);
		params.addm("Ry", 2160000.0);
		params.add("L", 6000000.0);
		params.add("dx", 30000.0);
		params.addm("ky", 2);
		params.addm("J", -1.0);
		params.add("Random seed", 0);
		params.add("Magnetization", 0.0);
		params.addm("range change", 0.01);
		params.addm("dt", 0.001);
		params.add("dt new");
		params.add("Time");
		params.add("Mean Phi");
		flags.add("Write 1D Config");
	}

	public void animate() {
		params.set("Time", ising.time());
		params.set("Mean Phi", ising.mean(ising.phi));
		phiGrid.registerData(Lp, Lp, ising.phi);
		etaVsTimeSim.setAutoScale(true);
		etaVsTimeLinear.setAutoScale(true);
		etaVsTimeLinearK.setAutoScale(true);
		etaVsTimeLC.setAutoScale(true);
		hSlice.setAutoScale(true);
		vSlice.setAutoScale(true);		
		
		double horizontalSlice = params.fget("Horizontal Slice");
		double verticalSlice = params.fget("Vertical Slice");
		
		hSlice.registerLines("Slice", ising.getHslice(horizontalSlice), Color.GREEN);
		hSlice.registerLines("phi0", new PointSet(0, 1, phi0) , Color.BLACK);
		vSlice.registerLines("Slice", ising.getVslice(verticalSlice), Color.BLUE);
		
		phiGrid.setDrawables(asList(
				Geom2D.line(0, horizontalSlice, 1, horizontalSlice, Color.GREEN),
				Geom2D.line(verticalSlice, 0, verticalSlice, 1, Color.BLUE)));
		
		if(accEtaValues){
			etaVsTimeSim.registerLines("acc", etaAcc, Color.BLACK);
			etaVsTimeSim.registerLines("acc2", etaAcc2, Color.BLUE);
			etaVsTimeSim.registerLines("acc3", etaAcc3, Color.RED);
			etaVsTimeSim.registerLines("acc4", etaAcc4, Color.ORANGE);
		
			etaVsTimeLinear.registerLines("accLT", etaLTAcc, Color.BLACK);
			etaVsTimeLinear.registerLines("acc2LT", etaLTAcc2, Color.BLUE);
			etaVsTimeLinear.registerLines("acc3LT", etaLTAcc3, Color.RED);
			etaVsTimeLinear.registerLines("acc4LT", etaLTAcc4, Color.ORANGE);

			etaVsTimeLinearK.registerLines("accLTk", etaLTkAcc, Color.BLACK);
			etaVsTimeLinearK.registerLines("acc2LTk", etaLTkAcc2, Color.BLUE);
			etaVsTimeLinearK.registerLines("acc3LTk", etaLTkAcc3, Color.RED);
			etaVsTimeLinearK.registerLines("acc4LTk", etaLTkAcc4, Color.ORANGE);

			etaVsTimeLC.registerLines("accLC", etaLCAcc, Color.BLACK);
			etaVsTimeLC.registerLines("acc2LC", etaLCAcc2, Color.BLUE);
			etaVsTimeLC.registerLines("acc3LC", etaLCAcc3, Color.RED);
			etaVsTimeLC.registerLines("acc4LC", etaLCAcc4, Color.ORANGE);
		}

		double [] eigenVect2D = new double [Lp*Lp];
		for (int i = 0; i < Lp; i++){
			for (int  j = 0; j < Lp; j++){
				eigenVect2D[i+Lp*j]=VV[j][i];
			}
		}

		if(flags.contains("Clear")){
			etaAcc.clear();
			etaAcc2.clear();
			etaAcc3.clear();
			etaAcc4.clear();

			etaLTAcc.clear();
			etaLTAcc2.clear();
			etaLTAcc3.clear();
			etaLTAcc4.clear();
			
			etaLTkAcc.clear();
			etaLTkAcc2.clear();
			etaLTkAcc3.clear();
			etaLTkAcc4.clear();
			
			etaLCAcc.clear();
			etaLCAcc2.clear();
			etaLCAcc3.clear();
			etaLCAcc4.clear();
			
			flags.clear();			
		}
		
	}

	public void clear() {
		etaAcc.clear();
		etaLTAcc.clear();
		etaAcc2.clear();
		etaLTAcc2.clear();
		etaAcc3.clear();
		etaLTAcc3.clear();
		etaAcc4.clear();
		etaLTAcc4.clear();
	}
 
	public void run() {
		clearFile = true;
		initialize();
		ky = params.iget("ky");
		double [] etaK = new double[Lp*Lp];
		double recordStep = .0001;	
		
		for (int i = 0; i < Lp*Lp; i++)
			etaLT[i] = ising.phi[i] - phi0[i%Lp];
		double [] etaLT_k = new double [Lp*Lp];
		etaLT_k = fft.calculate2DFT(etaLT);
		for (int i = ky*Lp; i < Lp*(ky+1); i++)
			etaLT_k_slice[i-ky*Lp] = etaLT_k[i]; 
		for (int i = 0; i < Lp*Lp; i++)
			eta[i] = ising.phi[i] - phi0[i%Lp];

		System.out.println("ky = " + ky);

		if(modifiedDynamics) findModMatrix();
		else findMatrix();

		diagonalize();
		findAndDiagonalizeLittleMatrix();
		Job.animate();

		while (true) {
			if (flags.contains("Write 1D Config")){
				write1Dconfig();
				flags.clear();
			}
			if(modifiedDynamics){
        		ising.simulate();
        		if(accEtaValues){
        			//rhs2D = simulateLinearWithModDynamics(etaLT);
        			//rhs = simulateLinearKwithModDynamics();
        		}
        	}else{
        		//System.out.println(ising.dt);
        		ising.readParams(params);
        		ising.simulateUnstable();
        		//System.out.println(ising.dt + " dt");
        		params.set("dt new", ising.dt);
       			if(accEtaValues){	
       				rhs2D = simulateLinear(etaLT);	
       				//rhs = simulateLinearKbar();
        			rhs = simulateLinearK();
        			etaLinearCombo();
       			}
        	}
			for (int i = 0; i < Lp*Lp; i++)
				eta[i] = ising.phi[i] - phi0[i%Lp];
			etaK = fft.calculate2DFT(eta);			
			if(accEtaValues){
				
				for (int i = 0; i < Lp*Lp; i++)
    				etaLT[i] += rhs2D[i];

				for (int i = 0; i < Lp; i++)
    				etaLT_k_slice[i] += rhs[i];

   				etaAcc.accum(ising.time(), Math.pow(etaK[ky*Lp],2));
				etaAcc2.accum(ising.time(),Math.pow(etaK[ky*Lp+1],2));
				etaAcc3.accum(ising.time(), Math.pow(etaK[ky*Lp+2],2));
				etaAcc4.accum(ising.time(), Math.pow(etaK[ky*Lp+3],2)); 
				
				double [] etaLTSF = new double[Lp*Lp];
				etaLTSF = fft.calculate2DFT(etaLT);
				etaLTAcc.accum(ising.time(),Math.pow(etaLTSF[ky*Lp],2));		
				etaLTAcc2.accum(ising.time(), Math.pow(etaLTSF[ky*Lp+1],2));
				etaLTAcc3.accum(ising.time(),Math.pow(etaLTSF[ky*Lp+2],2));
				etaLTAcc4.accum(ising.time(), Math.pow(etaLTSF[ky*Lp+3],2));
				
				etaLTkAcc.accum(ising.time(),Math.pow(etaLT_k_slice[0],2));		
				etaLTkAcc2.accum(ising.time(), Math.pow(etaLT_k_slice[1],2));
				etaLTkAcc3.accum(ising.time(),Math.pow(etaLT_k_slice[2],2));
				etaLTkAcc4.accum(ising.time(), Math.pow(etaLT_k_slice[3],2));
		
				etaLCAcc.accum(ising.time(),Math.pow(etaLC[0],2));		
				etaLCAcc2.accum(ising.time(), Math.pow(etaLC[1],2));
				etaLCAcc3.accum(ising.time(),Math.pow(etaLC[2],2));
				etaLCAcc4.accum(ising.time(), Math.pow(etaLC[3],2)); 					
			}
   		
    		Job.animate();
    		if(writeToFile){
    			if(ising.t >= recordStep){
    				recordSfDataToFile(etaK);
    				recordStep += .01;
    			}
    		}
		}	
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
			etaLC[i] = sum/(double)(Lp);
		}
	}
	
	/**
	* Evaluates the matrix, M, of the linear theory for the unstable 
	* (unmodified) Ising dynamics.
	*/
	private void findMatrix() {
		double kyValue = 2.0*Math.PI*ising.Ry*ky/ising.L;
		double kxValue;
		for (int i = 0; i < Lp; i++){
			for (int j = 0; j <= i; j++){
				M[i][j]=M[j][i]=-ising.T*f_k[(i-j+Lp)%Lp]/Lp;
			}
		}
		for (int i = 0; i < Lp; i++){
			if(i >= Lp/2)
				kxValue = 2.0*Math.PI*ising.Rx*(i-Lp)/ising.L;
			else
				kxValue = 2.0*Math.PI*ising.Rx*i/ising.L;
			if(ising.circleInteraction)
				M[i][i]-=ising.findVkCircle(Math.sqrt(kxValue*kxValue + kyValue*kyValue));
			else
				M[i][i]-=ising.findVkSquare(kxValue, kyValue);
		}
	}

	/**
	* Evaluates the matrix, M, of the linear theory for the stable 
	* (modified) Ising dynamics.  Not working at present.
	*/
	void findModMatrix(){
	}	
	
	private void findAndDiagonalizeLittleMatrix() {
		
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
		
	double [] simulateLinearK(){
		double [] linearTheoryGrowth = new double[Lp*Lp];
				
		for (int i = 0; i < Lp; i++){
			for (int j = 0; j < Lp; j++){
				linearTheoryGrowth [i] += ising.dt*M[i][j]*etaLT_k_slice[j];
			}
		}
 
		return linearTheoryGrowth;
	}

	double [] simulateLinearKbar(){
		double kyValue = 2.0*Math.PI*ising.Ry*ky/ising.L;
		double [] linearTheoryGrowth = new double[Lp*Lp];
		double kxValue;
		double [] f_bar = new double [Lp];
		f_bar = fft.backConvolve1D(f_k, etaLT_k_slice);
		convolve.registerLines("convolvution", new PointSet(1,1, f_bar), Color.green);
		for (int i = 0; i < Lp; i++){
			
			if(i >= Lp/2)
				kxValue = 2.0*Math.PI*ising.Rx*(i-Lp)/ising.L;
			else
				kxValue = 2.0*Math.PI*ising.Ry*i/ising.L;
				
			linearTheoryGrowth[i] = ising.dt*(-ising.findVkSquare(kxValue, kyValue)*etaLT_k_slice[i]);
			linearTheoryGrowth [i] -= ising.dt*(ising.T*f_bar[i]);
		}
		etaLTkAcc.accum(ising.time(),Math.pow(etaLT_k_slice[0],2));		
   		etaLTkAcc2.accum(ising.time(), Math.pow(etaLT_k_slice[1],2));
   		etaLTkAcc3.accum(ising.time(),Math.pow(etaLT_k_slice[2],2));
   		etaLTkAcc4.accum(ising.time(), Math.pow(etaLT_k_slice[3],2)); 
		return linearTheoryGrowth;
	}

	/**
	* Diagonalizes the matrix M.
	* Also some code that can be used to check for
	* orthonormality of eigenvectors. (Eigenvectors should
	* be orthonormal for Hermitian matrices.  Matrix for unmodified
	* dynamics is Hermitian, but matrix for modified dynamcis
	* is not.
	*/
	void diagonalize(){
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
		double [] normedEta = MathTools.normalize(etaLT_k_slice);
		double normCheck = 0;
		for (int i = 0; i < Lp; i ++){
			c[i] = MathTools.dot(normedEta, VV[i]);
			normCheck += c[i]*c[i];
		}
		
		double sum = 0;
		int testInt=1;
	    for(int i = 0; i < Lp; i++){
	    		sum += M[testInt][i]*VV[testInt][i];
	    } 
//	    double lambda = sum/VV[testInt][testInt];
//	    System.out.println("ev = " + eigenvalue[testInt] + " lambda = " + lambda);
	    
	}
	
	void findPhi0andPhi0_bar(){
		String fileName = "../../../research/javaData/configs1d/configOpt";
//		String fileName = "../../../research/javaData/configs1d/config";
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
				double kRx = 2*Math.PI*ising.Rx*k1/ising.L;
				double kRy = 0.0;
				if(ising.circleInteraction) return ising.findVkCircle(kRx*kRx + kRy*kRy);
				else return ising.findVkSquare(kRx, kRy);
			}
		});
		for(int i = 0; i < Lp; i++)
			f_x[i] = 1.0 / (1.0 - Math.pow(phi0[i],2));
		f_k = fft.calculate1DFT(f_x);
	}
		
	double [] simulateLinearMod_kSpace(int ky){
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
	
//	double [] simulateLinearWithModDynamics(double [] etaLinear){
//		double [] linearTheoryGrowth = new double [Lp*Lp];
//		double [] eta_bar = new double[Lp*Lp];
//		ising.convolveWithRange(etaLinear, eta_bar, ising.R);
//		for (int i = 0; i < Lp*Lp; i++)
//			linearTheoryGrowth[i] = +ising.dt*(eta_bar[i]*G[i%Lp]+ etaLinear[i]*F[i%Lp]);
//		return linearTheoryGrowth;
//	}
	
	double [] simulateLinear(double [] etaLinear){
		double [] linearTheoryGrowth = new double[Lp*Lp];
		double [] eta_bar = new double[Lp*Lp];
		double [] phi0_2D = new double[Lp*Lp];
		for(int i = 0; i < Lp*Lp; i++)
			phi0_2D[i] = phi0[i%Lp];
//		ising.convolveWithRange(eta, eta_bar, ising.R);
		eta_bar = fft.convolve2DwithFunction(etaLinear, new Function2D(){
			public double eval(double k1, double k2) {
				double kRx = 2*Math.PI*ising.Rx*k1/ising.L;
				double kRy = 2*Math.PI*ising.Ry*k2/ising.L;
				if(ising.circleInteraction) return ising.findVkCircle(kRx*kRx + kRy*kRy);
				else return ising.findVkSquare(kRx, kRy);
			}
		});

		for (int i = 0; i < Lp*Lp; i++)
			linearTheoryGrowth[i] = ising.dt*(-eta_bar[i]-ising.T*etaLinear[i]*f_x[i%Lp]);
		return linearTheoryGrowth;
	}
	
	public void recordSfDataToFile(double [] data){
		String file0 = "../../../research/javaData/stripeToClumpInvestigation/kySlice/sf0";
		String file1 = "../../../research/javaData/stripeToClumpInvestigation/kySlice/sf1";
		String file2 = "../../../research/javaData/stripeToClumpInvestigation/kySlice/sf2";
		String file3 = "../../../research/javaData/stripeToClumpInvestigation/kySlice/sf3";
		String file4 = "../../../research/javaData/stripeToClumpInvestigation/kySlice/sf4";
		String file5 = "../../../research/javaData/stripeToClumpInvestigation/kySlice/sf5";
		if (clearFile){
			FileUtil.deleteFile(file0);
			FileUtil.deleteFile(file1);
			FileUtil.deleteFile(file2);
			FileUtil.deleteFile(file3);
			FileUtil.deleteFile(file4);
			FileUtil.deleteFile(file5);
			clearFile = false;
		}
		FileUtil.printlnToFile(file0, ising.time(), data[ky*Lp]*data[ky*Lp]);
		FileUtil.printlnToFile(file1, ising.time(), data[ky*Lp+1]*data[ky*Lp+1]);
		FileUtil.printlnToFile(file2, ising.time(), data[ky*Lp+2]*data[ky*Lp+2]);
		FileUtil.printlnToFile(file3, ising.time(), data[ky*Lp+3]*data[ky*Lp+3]);
		FileUtil.printlnToFile(file4, ising.time(), data[ky*Lp+4]*data[ky*Lp+4]);
		FileUtil.printlnToFile(file5, ising.time(), data[ky*Lp+5]*data[ky*Lp+5]);
		System.out.println("data written for time = " + ising.time());
	}	
	
	public void readInputParams(String FileName){
		try {
			File inputFile = new File(FileName);
			DataInputStream dis = new DataInputStream(new FileInputStream(inputFile));
			double readData;
			
			readData = dis.readDouble();
			System.out.println(readData);
			params.set("J", readData);
			dis.readChar();				
			
			readData = dis.readDouble();
			System.out.println(readData);
			params.set("R", readData);
			dis.readChar();	
			
			readData = dis.readDouble();
			System.out.println(readData);
			params.set("L/R", readData);
			dis.readChar();
			
			readData = dis.readDouble();
			System.out.println(readData);
			params.set("R/dx", readData);
			dis.readChar();	
			
			System.out.println("input read");
		}catch(IOException ex){
			ex.printStackTrace();
		}
	}
	
	void initialize(){

		ising = new IsingField2Dopt(params);
//		Does a random config by default:
//		Have to read in from file separately
		if(params.sget("Init Conditions") == "Read From File")
			ising.readConfigFromFile();
//		readInputParams("../../../research/javaData/configs/inputParams");
		
		etaAcc = new Accumulator(ising.dt);
		etaAcc2 = new Accumulator(ising.dt);
		etaAcc3 = new Accumulator(ising.dt);
		etaAcc4 = new Accumulator(ising.dt);
		
		etaLTAcc = new Accumulator(ising.dt);		
		etaLTAcc2 = new Accumulator(ising.dt);
		etaLTAcc3 = new Accumulator(ising.dt);
		etaLTAcc4 = new Accumulator(ising.dt);
			

		etaLTkAcc = new Accumulator(ising.dt);		
		etaLTkAcc2 = new Accumulator(ising.dt);
		etaLTkAcc3 = new Accumulator(ising.dt);
		etaLTkAcc4 = new Accumulator(ising.dt);

		etaLCAcc = new Accumulator(ising.dt);		
		etaLCAcc2 = new Accumulator(ising.dt);
		etaLCAcc3 = new Accumulator(ising.dt);
		etaLCAcc4 = new Accumulator(ising.dt);
		
		this.Lp = ising.Lp;
		fft = new FourierTransformer(Lp);
		rhs = new double[Lp];
		rhs2D = new double[Lp*Lp];
		eta = new double[Lp*Lp];
		etaLT = new double[Lp*Lp];
		//etaLT_k = new double[Lp*Lp];
		etaLT_k_slice = new double [Lp];
		phi0 = new double [Lp];
		phi0_bar = new double [Lp];
		g_k = new double[Lp];
		f_k = new double[Lp];
		f_x = new double [Lp];
		F = new double[Lp];
		G = new double[Lp];
		M = new double [Lp][Lp];
		c = new double[Lp];
		etaLC = new double[Lp];
		eigenvalue = new double[Lp];
		VV = new double [Lp][Lp];
		findPhi0andPhi0_bar();

		if(params.sget("Init Conditions") == "Read in 1D solution"){
			System.out.println("Read in 1D solution");
			for(int i = 0; i < Lp*Lp; i++)
				ising.phi[i]=phi0[i%Lp] + (random.nextGaussian())/1000000;
		}
	}	

	private void write1Dconfig(){
		String configFileName = "../../../research/javaData/configs1d/configOpt";
		FileUtil.deleteFile(configFileName);
		FileUtil.writeConfigToFile(configFileName, ising.Lp, ising.phi);
	}
	
}


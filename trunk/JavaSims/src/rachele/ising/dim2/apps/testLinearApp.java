package rachele.ising.dim2.apps;

import java.io.DataInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.awt.Color;
import rachele.ising.dim2.FindCoefficients;
import rachele.util.FourierTransformer;
import rachele.util.MathTools;
import rachele.ising.dim2.IsingField2D;
import rachele.util.FileUtil;
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


public class testLinearApp extends Simulation{
    IsingField2D ising;
    public FindCoefficients coeff;
    public FourierTransformer fft;
    public double [] rhs, rhs2D, etaLT, eta, etaLT_k_slice, phi0_bar, phi0; //right hand side
    public double [] g_k, f_x, f_k, f, g, c, eigenvalue, etaLC;
    public double [][] M, VV;
    public int ky;
    public double kRChunk; //=2piR/L

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
    Grid evGrid = new Grid("eigenvectors");
    Plot evPlot = new Plot("eigenvalues");
    Plot etaVsTimeSim = new Plot("eta v t");
    Plot etaVsTimeLinearK = new Plot("eta v t LT k space");
    Plot etaVsTimeLinear = new Plot("eta v t LT");
    Plot etaVsTimeLC = new Plot("eta v t LC");
    
    Plot fkPlot = new Plot ("f(k)");
    Plot etakSimPlot = new Plot("eta k sim");
    Plot convolve = new Plot("convolve");
    //Plot  = new Plot("g(k)");
    
    
	public static void main(String[] args) {
		new Control(new testLinearApp(), "Ising Linear Test");
	}
	public void load(Control c) {
		c.frameTogether("accs", etaVsTimeSim, etaVsTimeLinear, etaVsTimeLinearK, etaVsTimeLC);
		//c.frameTogether("checks", fkPlot, etakPlot, etakSimPlot, convolve);
		c.frameTogether("eigens",evGrid,evPlot);
		params.addm("Zoom", new ChoiceValue("Yes", "No"));
		params.addm("Interaction", new ChoiceValue("Square", "Circle"));
		params.addm("Dynamics?", new ChoiceValue("Langevin No M Convervation"));
		params.add("Init Conditions", new ChoiceValue("Read 1D Soln", "Read From File","Random Gaussian"));
		params.addm("Approx", new ChoiceValue("Slow", "HalfStep", "TimeAdjust", "Phi4","Phi4HalfStep"));
		params.addm("Noise", new DoubleValue(0.0, 0.0, 1.0).withSlider());
		params.addm("Horizontal Slice", new DoubleValue(0.5, 0, 0.9999).withSlider());
		params.addm("Vertical Slice", new DoubleValue(0.5, 0, 0.9999).withSlider());
		params.addm("kR", new DoubleValue(5.135622302, 0.0, 6.0).withSlider());
		params.addm("T", 0.04);
		params.addm("H", 0.80);
		params.addm("dT", 0.001);
		params.addm("tolerance", 0.0001);
		params.addm("dt", 0.01);
		params.addm("J", -1.0);
		params.addm("R", 2000000.0);
		params.addm("Random seed", 0);
		params.add("L/R", 3.0);
		params.add("R/dx", 43.0);
		params.add("kR bin-width", 0.1);
		params.add("Magnetization", 0.0);
		params.addm("ky", 2);
		params.add("Time");
		params.add("Mean Phi");
		params.add("Lp");
		flags.add("Clear");
	}

	public void animate() {
		params.set("Time", ising.time());
		params.set("Mean Phi", ising.mean(ising.phi));
		params.set("Lp", ising.Lp);
		phiGrid.registerData(Lp, Lp, ising.phi);
		etaVsTimeSim.setAutoScale(true);
		etaVsTimeLinear.setAutoScale(true);
		etaVsTimeLinearK.setAutoScale(true);
		etaVsTimeLC.setAutoScale(true);
		
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
		
		//etakPlot.registerLines("eta(k) check", new PointSet(1,1,etaLT_k_slice), Color.BLACK);
		fkPlot.registerLines("f(k) check", new PointSet(1,1,f_k), Color.BLUE);
		fkPlot.registerLines("f(x) check", new PointSet(1,1,f_x), Color.YELLOW);
		double [] tester = new double [Lp];
		tester = fft.calculate1DBackFT(f_k);
		fkPlot.registerLines("f(x) check 2", new PointSet(1,1,tester), Color.BLACK);
		double [] etaSF = new double[Lp*Lp];
		etaSF = fft.calculate2DFT(eta);
		double [] etaSF_slice = new double [Lp];
		for (int i = ky*Lp; i < Lp*(ky+1); i++)
			etaSF_slice[i-ky*Lp] = etaSF[i]*etaSF[i];
		etakSimPlot.registerLines("eta(k) sim check", new PointSet(1,1,etaSF_slice), Color.RED);
		
		double [] eigenVect2D = new double [Lp*Lp];
		for (int i = 0; i < Lp; i++){
			for (int  j = 0; j < Lp; j++){
				eigenVect2D[i+Lp*j]=VV[j][i];
			}
		}
		evGrid.registerData(Lp,Lp,eigenVect2D);
		evPlot.registerLines("eigenvalues", new PointSet(1,1,eigenvalue), Color.BLACK);
		
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
		initialize();
		ky = params.iget("ky");
		for (int i = 0; i < Lp*Lp; i++)
			etaLT[i] = ising.phi[i] - phi0[i%Lp];
		double [] etaLT_k = new double [Lp*Lp];
		etaLT_k = fft.calculate2DFT(etaLT);
		for (int i = ky*Lp; i < Lp*(ky+1); i++)
			etaLT_k_slice[i-ky*Lp] = etaLT_k[i]; 

		for (int i = 0; i < Lp*Lp; i++)
			eta[i] = ising.phi[i] - phi0[i%Lp];
		boolean modifiedDynamics = false;
//		boolean realSpace = false;
		if(modifiedDynamics){
			findModMatrix();
			diagonalize();
		}else{
			findMatrix();
			diagonalize();
		}
//		int trialInt = 10;
//		for (int i = 0; i < Lp; i++)
//			etaLT_k_slice[i] = VV[trialInt][i];
//		System.out.println("ev " + trialInt + " = " + eigenvalue[trialInt]);
		Job.animate();
		while (true) {
        	if(modifiedDynamics){
        		rhs2D = simulateLinearWithModDynamics(etaLT);
        		//simulateLinearMod_kSpace();
        		//rhs = simulateLinearMod_kSpace(ky);
        		ising.simulate();
        		for (int i = 0; i < Lp*Lp; i++){
        			eta[i] = ising.phi[i] - phi0[i%Lp];
        			etaLT[i] += rhs2D[i];
        		}
        	}else{
        		ising.simulateUnstable();
    			for (int i = 0; i < Lp*Lp; i++)
    				eta[i] = ising.phi[i] - phi0[i%Lp];
    			double [] etaSF = new double[Lp*Lp];
    			etaSF = fft.calculate2DFT(eta);
    			etaAcc.accum(ising.time(), Math.pow(etaSF[ky*Lp],2));
    			etaAcc2.accum(ising.time(),Math.pow(etaSF[ky*Lp+1],2));
    			etaAcc3.accum(ising.time(), Math.pow(etaSF[ky*Lp+2],2));
    			etaAcc4.accum(ising.time(), Math.pow(etaSF[ky*Lp+3],2)); 
      		
//        		if(realSpace){
        			rhs2D = simulateLinear(etaLT);	
        			for (int i = 0; i < Lp*Lp; i++)
        				etaLT[i] += rhs2D[i];
        			double [] etaLTSF = new double[Lp*Lp];
        			etaLTSF = fft.calculate2DFT(etaLT);
//        			for (int i = ky*Lp; i < Lp*(ky+1); i++)
//        				etaLT_k_slice[i-ky*Lp] = etaLTSF[i]; 
        			etaLTAcc.accum(ising.time(),Math.pow(etaLTSF[ky*Lp],2));		
        			etaLTAcc2.accum(ising.time(), Math.pow(etaLTSF[ky*Lp+1],2));
        			etaLTAcc3.accum(ising.time(),Math.pow(etaLTSF[ky*Lp+2],2));
        			etaLTAcc4.accum(ising.time(), Math.pow(etaLTSF[ky*Lp+3],2));
        			//etakPlot.registerLines("eta(k) check", new PointSet(1,1,etaLT_k_slice), Color.BLACK);
//        		}else{
        			//rhs = simulateLinearKbar();
        				rhs = simulateLinearK();
        				for (int i = 0; i < Lp; i++)
        					etaLT_k_slice[i] += rhs[i];
        				etakPlot.registerLines("eta(k) check", new PointSet(1,1,etaLT_k_slice), Color.BLACK);
        				etaLinearCombo();
        		//}
        	}
    		Job.animate();
		}	
	}
	
	public void etaLinearCombo(){
		
		for(int i = 0; i < Lp; i++){
			double sum = 0;
			for(int j = 0; j < Lp; j++){
				sum += c[j]*VV[j][i]*Math.exp(ising.time()*eigenvalue[j]);
			}
			etaLC[i] = sum;
		}
		etaLCAcc.accum(ising.time(),Math.pow(etaLC[0],2));		
		etaLCAcc2.accum(ising.time(), Math.pow(etaLC[1],2));
		etaLCAcc3.accum(ising.time(),Math.pow(etaLC[2],2));
		etaLCAcc4.accum(ising.time(), Math.pow(etaLC[3],2)); 		
	}
	
	private void findMatrix() {
		double kyValue = 2.0*Math.PI*ising.R*ky/ising.L;
		double kxValue;
		for (int i = 0; i < Lp; i++){
			for (int j = 0; j <= i; j++){
				//M[i][j]=-ising.T*f_k[(i-j+Lp)%Lp]/Lp;
				M[i][j]=M[j][i]=-ising.T*f_k[(i-j+Lp)%Lp]/Lp;
			}
		}
		for (int i = 0; i < Lp; i++){
			if(i >= Lp/2)
				kxValue = 2.0*Math.PI*ising.R*(i-Lp)/ising.L;
			else
				kxValue = 2.0*Math.PI*ising.R*i/ising.L;
			M[i][i]-=ising.findVkSquare(kxValue, kyValue);
		}
	}

	
	double [] simulateLinearK(){
		double [] linearTheoryGrowth = new double[Lp*Lp];
				
		for (int i = 0; i < Lp; i++){
			for (int j = 0; j < Lp; j++){
				linearTheoryGrowth [i] += ising.dt*M[i][j]*etaLT_k_slice[j];
				//linearTheoryGrowth [i] += ising.dt*M[i][j]*VV[0][j];
			}
		}
			etaLTkAcc.accum(ising.time(),Math.pow(etaLT_k_slice[0],2));		
   			etaLTkAcc2.accum(ising.time(), Math.pow(etaLT_k_slice[1],2));
   			etaLTkAcc3.accum(ising.time(),Math.pow(etaLT_k_slice[2],2));
   			etaLTkAcc4.accum(ising.time(), Math.pow(etaLT_k_slice[3],2)); 
		return linearTheoryGrowth;
	}

	double [] simulateLinearKbar(){
		double kyValue = 2.0*Math.PI*ising.R*ky/ising.L;
		double [] linearTheoryGrowth = new double[Lp*Lp];
		double kxValue;
		double [] f_bar = new double [Lp];
		f_bar = fft.backConvolve1D(f_k, etaLT_k_slice);
		convolve.registerLines("convolvution", new PointSet(1,1, f_bar), Color.green);
		for (int i = 0; i < Lp; i++){
			
			if(i >= Lp/2)
				kxValue = 2.0*Math.PI*ising.R*(i-Lp)/ising.L;
			else
				kxValue = 2.0*Math.PI*ising.R*i/ising.L;
				
			linearTheoryGrowth[i] = ising.dt*(-ising.findVkSquare(kxValue, kyValue)*etaLT_k_slice[i]);
			linearTheoryGrowth [i] -= ising.dt*(ising.T*f_bar[i]);
		}
			etaLTkAcc.accum(ising.time(),Math.pow(etaLT_k_slice[0],2));		
   			etaLTkAcc2.accum(ising.time(), Math.pow(etaLT_k_slice[1],2));
   			etaLTkAcc3.accum(ising.time(),Math.pow(etaLT_k_slice[2],2));
   			etaLTkAcc4.accum(ising.time(), Math.pow(etaLT_k_slice[3],2)); 
		return linearTheoryGrowth;
	}

	
	void diagonalize(){
		Matrix matrix = new Matrix(M);
		EigenvalueDecomposition Eig;
		Eig = matrix.eig();
		eigenvalue = new double [Lp];
		eigenvalue = Eig.getRealEigenvalues();
		
		SingularValueDecomposition Svd = matrix.svd();
		double [] singularValue = new double [Lp];
		singularValue = Svd.getSingularValues();
		
			
//		for (int i = 0; i < Lp; i ++){
//			System.out.println("eigenvalue " + i + " = " + eigenvalue[i]);
//			System.out.println("singular value " + i + " = " + singularValue[i]);
//		}
		System.out.println("eigenvalue max ="  +  DoubleArray.max(eigenvalue));
		System.out.println("eigenvalue min ="  + DoubleArray.min(eigenvalue));
		System.out.println("Svalue max ="  +  DoubleArray.max(singularValue));
		System.out.println("Svalue min ="  + DoubleArray.min(singularValue));
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
			//System.out.println("norm check = " + i + " = " + c[i]);
		}
		
		double sum = 0;
		int testInt=1;
	    for(int i = 0; i < Lp; i++){
	    		sum += M[testInt][i]*VV[testInt][i];
	    } 
	    double lambda = sum/VV[testInt][testInt];
	    System.out.println("ev = " + eigenvalue[testInt] + " lambda = " + lambda);
	    
		//System.out.println("norm check = " + normCheck);
		
		//		Matrix V = Eig.getV();
//		Matrix D = Eig.getD();
//		double[] finalDvector = D.transpose().getArray()[Lp - 1];
//		double[] finalEigenvector = V.transpose().getArray()[Lp - 1];
//		for (int i = 0; i < Lp; i ++){
//			System.out.println(i + " " + finalEigenvector[i]);
//		}
//		System.out.println("eigenvector ? " + finalDvector[Lp - 1]);
	}
	
	void findPhi0andPhi0_bar(){
		String fileName = "../../../research/javaData/configs1d/config";
		phi0 = FileUtil.readConfigFromFile(fileName, Lp);
		phi0_bar = fft.convolve1DwithFunction(phi0, new Function1D(){
			public double eval(double k1) {
				double kRx = 2*Math.PI*ising.R*k1/ising.L;
				return ising.findVkSquare(kRx, 0);
			}
		});
		for(int i = 0; i < Lp; i++)
			f_x[i] = 1.0 / (1.0 - Math.pow(phi0[i],2));
		f_k = fft.calculate1DFT(f_x);
	}
		
	void findModMatrix(){
		for (int i = 0; i < Lp; i++){
			double alpha_x = 1.0-Math.pow(phi0[i],2);
			g[i] = -alpha_x*alpha_x;
			f[i] = alpha_x*(ising.T*(4*phi0[i]*scikit.numerics.Math2.atanh(phi0[i])-1.0)+4.0*phi0[i]*(phi0_bar[i]-ising.H));				
		}
		g_k = fft.calculate1DFT(g);
		f_k = fft.calculate1DFT(f);

		//for (int ky = -Lp/2; ky < Lp/2; ky++){
		double kyValue = 2.0*Math.PI*ising.R*ky/ising.L;
		for(int kx = -Lp/2; kx < Lp/2; kx++){	
			//int i = Lp*((ky+Lp)%Lp) + (kx+Lp)%Lp;
			//double sum = 0;				
			for (int kxx = -Lp/2; kxx < Lp/2; kxx++){
				double kxxValue = 2.0*Math.PI*ising.R*kxx/ising.L;
				//int eta_k_int = Lp*((ky+Lp)%Lp) + (kxx+Lp)%Lp;
				M[(kx+Lp)%Lp][(kxx+Lp)%Lp] += (g_k [(kxx-kx+Lp)%Lp] * ising.findVkSquare(kxxValue,kyValue) + f_k[(kxx-kx+Lp)%Lp]);				
			}
		}
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
	
	double [] simulateLinearWithModDynamics(double [] etaLinear){
		double [] linearTheoryGrowth = new double [Lp*Lp];
		double [] eta_bar = new double[Lp*Lp];
		ising.convolveWithRange(etaLinear, eta_bar, ising.R);
		for (int i = 0; i < Lp*Lp; i++)
			linearTheoryGrowth[i] = +eta_bar[i]*g[i%Lp] + etaLinear[i]*f[i%Lp];
		return linearTheoryGrowth;
	}
	
	double [] simulateLinear(double [] etaLinear){
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
				return ising.findVkSquare(kRx, kRy);
			}
		});

		for (int i = 0; i < Lp*Lp; i++)
			linearTheoryGrowth[i] = ising.dt*(-eta_bar[i]-ising.T*etaLinear[i]*f_x[i%Lp]);
		return linearTheoryGrowth;
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
		if(params.sget("Init Conditions") == "Read From File")
			readInputParams("../../../research/javaData/configs/inputParams");
		ising = new IsingField2D(params);
		
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
		coeff = new FindCoefficients(Lp);
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
		f = new double[Lp];
		g = new double[Lp];
		M = new double [Lp][Lp];
		c = new double[Lp];
		etaLC = new double[Lp];
		eigenvalue = new double[Lp];
		VV = new double [Lp][Lp];
		findPhi0andPhi0_bar();
	}	

}

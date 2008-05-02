package rachele.ising.dim2.apps;

import java.io.DataInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.awt.Color;
import rachele.ising.dim2.FindCoefficients;
import rachele.util.FourierTransformer;
import rachele.ising.dim2.IsingField2D;
import rachele.util.FileUtil;
import scikit.graphics.dim2.Grid;
import scikit.graphics.dim2.Plot;
import scikit.jobs.Control;
import scikit.jobs.Job;
import scikit.jobs.Simulation;
import scikit.jobs.params.ChoiceValue;
import scikit.jobs.params.DoubleValue;
import scikit.numerics.Jama.EigenvalueDecomposition;
import scikit.numerics.Jama.Matrix;
import scikit.numerics.fn.Function1D;
import scikit.dataset.Accumulator;

public class testLinearApp extends Simulation{
    IsingField2D ising;
    public FindCoefficients coeff;
    public FourierTransformer fft;
    public double [] rhs, rhs2D, etaLT, eta, etaLT_k_slice, phi0_bar, phi0; //right hand side
    public double [] g_k, f_k, f, g;
    public double [][] M;
    public double kRChunk; //=2piR/L
    public Accumulator etaAcc;
    public Accumulator etaLTAcc; 
    public Accumulator etaAcc2;
    public Accumulator etaLTAcc2;
    public Accumulator etaAcc3;
    public Accumulator etaLTAcc3;
    public Accumulator etaAcc4;
    public Accumulator etaLTAcc4;
    public int Lp;
    Grid etaDot = new Grid("ising delta phi");
    Grid eta_k_LTGrid = new Grid("rhs");
    Grid phiGrid = new Grid("phi field");
    Grid etaDotSF = new Grid("sf phi field");
    Grid etaDotLinearSF = new Grid("sf linear theory");
    Grid phi0Grid = new Grid("phi_0");
    Grid phi0SliceGrid = new Grid("phi_0 slice");
    Plot etaVt1 = new Plot("eta v t");
    Plot etaVt2 = new Plot("eta v t LT");
    Plot fPlot = new Plot ("f(k)");
    Plot gPlot = new Plot("g(k)");
    
    
	public static void main(String[] args) {
		new Control(new testLinearApp(), "Ising Linear Test");
	}
	public void load(Control c) {
		c.frameTogether("accs", etaVt1, etaVt2);
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
		etaVt1.registerLines("acc", etaAcc, Color.BLACK);
		etaVt1.registerLines("acc2", etaAcc2, Color.BLUE);
		etaVt1.registerLines("acc3", etaAcc3, Color.RED);
		etaVt1.registerLines("acc4", etaAcc4, Color.ORANGE);
		
		etaVt2.registerLines("accLT", etaLTAcc, Color.BLACK);
		etaVt2.registerLines("acc2LT", etaLTAcc2, Color.BLUE);
		etaVt2.registerLines("acc3LT", etaLTAcc3, Color.RED);
		etaVt2.registerLines("acc4LT", etaLTAcc4, Color.ORANGE);
				
		if(flags.contains("Clear")){
			etaAcc.clear();
			etaLTAcc.clear();
			etaAcc2.clear();
			etaLTAcc2.clear();
			etaAcc3.clear();
			etaLTAcc3.clear();
			etaAcc4.clear();
			etaLTAcc4.clear();
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
		int ky = 4;
		for (int i = 0; i < Lp*Lp; i++)
			etaLT[i] = ising.phi[i] - phi0[i%Lp];
		double [] temp = new double [Lp*Lp];
		temp = fft.calculate2DFT(etaLT);
		for (int i = ky*Lp; i < Lp*(ky+1); i++)
			etaLT_k_slice[i-ky*Lp] = temp[i]; 
		boolean modifiedDynamics = false;
		boolean realSpace = true;
		if(modifiedDynamics){
			findModMatrix(ky);
			//diagonalize();
		}else{
			findMatrix(ky);
		}
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
        		if(realSpace){
        			rhs2D = simulateLinear(etaLT);
        			ising.simulateUnstable();
        			for (int i = 0; i < Lp*Lp; i++){
        				eta[i] = ising.phi[i] - phi0[i%Lp];
        				etaLT[i] += rhs2D[i];
        			}
        			double [] etaLTSF = new double[Lp*Lp];
        			etaLTSF = fft.calculate2DFT(etaLT);
        			etaLTAcc.accum(ising.time(),Math.pow(etaLTSF[ky*Lp],2));		
        			etaLTAcc2.accum(ising.time(), Math.pow(etaLTSF[ky*Lp+1],2));
        			etaLTAcc3.accum(ising.time(),Math.pow(etaLTSF[ky*Lp+2],2));
        			etaLTAcc4.accum(ising.time(), Math.pow(etaLTSF[ky*Lp+3],2));
        			double [] etaSF = new double[Lp*Lp];
        			etaSF = fft.calculate2DFT(eta);
        			etaAcc.accum(ising.time(), Math.pow(etaSF[ky*Lp],2));
        			etaAcc2.accum(ising.time(),Math.pow(etaSF[ky*Lp+1],2));
        			etaAcc3.accum(ising.time(), Math.pow(etaSF[ky*Lp+2],2));
        			etaAcc4.accum(ising.time(), Math.pow(etaSF[ky*Lp+3],2));
        		
        		}else{
        			ising.simulateUnstable();
        			for (int i = 0; i < Lp*Lp; i++){
        				eta[i] = ising.phi[i] - phi0[i%Lp];
        			}
        			
        			rhs = simulateLinearK(ky, etaLT_k_slice);
        			for (int i = 0; i < Lp; i++)
        				etaLT_k_slice[i] = rhs[i];
       				etaLTAcc.accum(ising.time(),Math.pow(etaLT_k_slice[0],2));		
           			etaLTAcc2.accum(ising.time(), Math.pow(etaLT_k_slice[1],2));
           			etaLTAcc3.accum(ising.time(),Math.pow(etaLT_k_slice[2],2));
           			etaLTAcc4.accum(ising.time(), Math.pow(etaLT_k_slice[3],2)); 
        			double [] etaSF = new double[Lp*Lp];
        			etaSF = fft.calculate2DFT(eta);
        			etaAcc.accum(ising.time(), Math.pow(etaSF[ky*Lp],2));
        			etaAcc2.accum(ising.time(),Math.pow(etaSF[ky*Lp+1],2));
        			etaAcc3.accum(ising.time(), Math.pow(etaSF[ky*Lp+2],2));
        			etaAcc4.accum(ising.time(), Math.pow(etaSF[ky*Lp+3],2));           			
        			
        		}
        	}
    		Job.animate();
		}	
	}
	
	private void findMatrix(int ky) {
		//double kyValue = 2.0*Math.PI*ising.R*ky/ising.L;
		double [] f_x = new double [Lp];
		for(int i = 0; i < Lp; i++)
			f_x[i] = 1.0 / (1.0 - Math.pow(phi0[i],2));
		f_k = fft.calculate1DFT(f_x);
		for(int kx = 0; kx < Lp; kx++){
			for(int kxx = 0; kxx < Lp; kxx++)
				M[kx][kxx] = -ising.T*f_k[(kx - kxx + Lp)%Lp];
		}
//		for(int kx = -Lp/2; kx < Lp/2; kx++){
//			double kxValue = 2.0*Math.PI*ising.R*kx/ising.L;
//			M[(kx+Lp)/Lp][(kx+Lp)/Lp] -= ising.findVkSquare(kxValue, kyValue);
//		}
	}

	double [] simulateLinearK(int ky, double [] etaK){
		double kyValue = 2.0*Math.PI*ising.R*ky/ising.L;
		double [] linearTheoryGrowth = new double[Lp*Lp];
		double kxValue;
		for (int i = 0; i < Lp; i++){
			
			if(i >= Lp/2)
				kxValue = 2.0*Math.PI*ising.R*(i-Lp)/ising.L;
			else
				kxValue = 2.0*Math.PI*ising.R*i/ising.L;
				
			linearTheoryGrowth[i] = ising.dt*(-ising.findVkSquare(kxValue, kyValue)*eta[i]);
			double sum = 0.0;
			for (int j = 0; j < Lp; j++)	
				sum += ising.dt*(M[i][j]*etaK[j]);
			linearTheoryGrowth [i] += sum/Lp;
		}
		return linearTheoryGrowth;
	}
	
	void diagonalize(){
		Matrix matrix = new Matrix(M);
		EigenvalueDecomposition Eig;
		Eig = matrix.eig();
		double [] eigenvalue = new double [Lp];
		eigenvalue = Eig.getRealEigenvalues();
		for (int i = 0; i < Lp; i ++)
			System.out.println("eigenvalue " + i + " = " + eigenvalue[i]);
		Matrix V = Eig.getV();
		Matrix D = Eig.getD();
		double[] finalDvector = D.transpose().getArray()[Lp - 1];
		double[] finalEigenvector = V.transpose().getArray()[Lp - 1];
		for (int i = 0; i < Lp; i ++){
			System.out.println(i + " " + finalEigenvector[i]);
		}
		System.out.println("eigenvector ? " + finalDvector[Lp - 1]);
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
	}
		
	void findModMatrix(int ky){
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
		ising.convolveWithRange(etaLinear, eta_bar, ising.R);
		for (int i = 0; i < Lp*Lp; i++)
			linearTheoryGrowth[i] = ising.dt*(+eta_bar[i]-ising.T*etaLinear[i]/(1.0-Math.pow(phi0[i%Lp],2)));
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
		etaLTAcc = new Accumulator(ising.dt);
		etaAcc2 = new Accumulator(ising.dt);
		etaLTAcc2 = new Accumulator(ising.dt);
		etaAcc3 = new Accumulator(ising.dt);
		etaLTAcc3 = new Accumulator(ising.dt);
		etaAcc4 = new Accumulator(ising.dt);
		etaLTAcc4 = new Accumulator(ising.dt);
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
		f = new double[Lp];
		g = new double[Lp];
		M = new double [Lp][Lp];
		findPhi0andPhi0_bar();
	}	

}

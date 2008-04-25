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
import scikit.numerics.fn.Function1D;
import scikit.dataset.Accumulator;
import scikit.dataset.PointSet;
import scikit.util.*;

public class testLinearApp extends Simulation{
    IsingField2D ising;
    public FindCoefficients coeff;
    public FourierTransformer fft;
    public double [] rhs, etaLT, eta, etaLT_k, phi0_bar, phi0; //right hand side
    public double [] g_k, f_k, f, g;
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
    //Grid etaDotLinear = new Grid("linear theory");
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
	c.frameTogether("Grids", etaDot, etaDotSF, phiGrid, eta_k_LTGrid, etaDotLinearSF);
	c.frameTogether("accs", etaVt1, etaVt2);
	c.frameTogether("f and g",fPlot, gPlot);
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
	params.addm("dt", 1.0);
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
		etaDot.registerData(Lp, Lp, ising.delPhi);
		eta_k_LTGrid.registerData(Lp, Lp, rhs);
		phiGrid.registerData(Lp, Lp, ising.phi);
		double [] delPhiFT = fft.calculateSF2D(ising.delPhi, false, false);
		double [] rhsFT = fft.calculateSF2D(rhs, false, false);
		double factor = delPhiFT [2*Lp]/DoubleArray.max(delPhiFT);
		double factorRHS = rhsFT[2*Lp]/DoubleArray.max(rhs);
		for (int i = 0; i < Lp; i++){
			delPhiFT[i] *= factor;
			rhsFT[i] *= factorRHS;
		}		

		fft.center(delPhiFT);
		fft.center(rhsFT);
		etaDotSF.registerData(ising.Lp, ising.Lp, delPhiFT);

		
		fft.center(delPhiFT);
		fft.center(rhsFT);
		double [] delPhiMidSlice = new double [Lp];
		double [] rhsMidSlice = new double [Lp];
		double [] delPhiSlice = new double [Lp];
		double [] rhsSlice = new double [Lp];
		for(int i = 0; i < Lp; i++){
			delPhiMidSlice[i] = delPhiFT[i];
			rhsMidSlice[i] = rhsFT[i];
			delPhiSlice[i] = delPhiFT[i*Lp+2];
			rhsSlice[i] = rhsFT[i*Lp+2];
		}

		double [] etaSF = new double[Lp*Lp];
		//double [] etaLTSF = new double[Lp*Lp];
		etaSF = fft.calculate2DFT(eta);
		//etaDotLinearSF.registerData(ising.Lp, ising.Lp, etaSF);
		//etaLTSF = fft.calculateSF2D(etaLT, false, false);
		etaLTAcc.accum(ising.time(),Math.pow(etaLT_k[2*Lp],2));
		etaAcc.accum(ising.time(), Math.pow(etaSF[2*Lp],2));
		etaLTAcc2.accum(ising.time(), Math.pow(etaLT_k[2*Lp+1],2));
		etaAcc2.accum(ising.time(),Math.pow(etaSF[2*Lp + 1],2));
		etaLTAcc3.accum(ising.time(),Math.pow(etaLT_k[2*Lp+2],2));
		etaAcc3.accum(ising.time(), Math.pow(etaSF[2*Lp + 2],2));
		etaLTAcc4.accum(ising.time(), Math.pow(etaLT_k[2*Lp+3],2));
		etaAcc4.accum(ising.time(), Math.pow(etaSF[2*Lp + 3],2));
		etaVt1.registerLines("acc", etaAcc, Color.BLACK);
		etaVt2.registerLines("acc", etaLTAcc, Color.BLACK);
		etaVt1.registerLines("acc2", etaAcc2, Color.BLUE);
		etaVt2.registerLines("acc2", etaLTAcc2, Color.BLUE);
		etaVt1.registerLines("acc3", etaAcc3, Color.RED);
		etaVt2.registerLines("acc3", etaLTAcc3, Color.RED);
		etaVt1.registerLines("acc4", etaAcc4, Color.ORANGE);
		etaVt2.registerLines("acc4", etaLTAcc4, Color.ORANGE);
		//fPlot.registerLines("f(k)", new PointSet(1,1,f_k), Color.BLACK);
		//fPlot.registerLines("g(k)", new PointSet(1,1,g_k), Color.BLUE);
		//fPlot.registerLines("phi0_k", new PointSet(1,1,fft.calculate1DFT(phi0)),Color.RED);
		//gPlot.registerLines("f(x)", new PointSet(1,1,f), Color.BLACK);
		//gPlot.registerLines("g(x)", new PointSet(1,1,g), Color.BLUE);
		//gPlot.registerLines("phi0", new PointSet(1,1,phi0), Color.RED);
		double [] etaSFslice = new double [Lp];
		double [] etaLTslice = new double[Lp];
		for (int i = 0; i < Lp; i++){
			etaSFslice[i] = etaSF[2*Lp+i];
			etaLTslice[i] = etaLT_k[2*Lp+i];
		}
		fPlot.registerLines("etaSlice", new PointSet(1,1,etaSFslice), Color.YELLOW);
		gPlot.registerLines("etaSlice2", new PointSet(1,1,etaLTslice), Color.GREEN);
		
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
		for (int i = 0; i < Lp*Lp; i++){
			etaLT[i] = ising.phi[i] - phi0[i%Lp];
		}
		etaLT_k = fft.calculate2DFT(etaLT);
		findGF();
        while (true) {
    		//rhs = simulateLinearWithModDynamics(etaLT);
    		//rhs = simulateLinear(phi0, etaLT);
        	//simulateLinearMod_kSpace();
       		rhs = simulateLinearMod_kSpace();
        	ising.simulate();
    		for (int i = 0; i < Lp*Lp; i++){
    			//etaLT[i] += rhs[i];
    			etaLT_k[i] += rhs[i];
    			eta[i] = ising.phi[i] - phi0[i%Lp];
    		}
    		Job.animate();
		}	
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
		
	void findGF(){
		for (int i = 0; i < Lp; i++){
			double alpha_x = 1.0-Math.pow(phi0[i],2);
			g[i] = -alpha_x*alpha_x;
			f[i] = alpha_x*(ising.T*(4*phi0[i]*scikit.numerics.Math2.atanh(phi0[i])-1.0)+4.0*phi0[i]*(phi0_bar[i]-ising.H));				
		}
		g_k = fft.calculate1DFT(g);
		f_k = fft.calculate1DFT(f);
	}
	
	double [] simulateLinearMod_kSpace(){
		double [] change = new double [Lp*Lp];

//		for (int ky = -Lp/2; ky < Lp/2; ky++){
		for (int ky = 0; ky < Lp; ky++){
		double kyValue = 2.0*Math.PI*ising.R*ky/ising.L;
			for (int kx = 0; kx < Lp; kx++){
				int i = Lp*((ky+Lp)%Lp) + (kx+Lp)%Lp;
				double sum = 0;				
				for (int kxx = 0; kxx < Lp; kxx++){
					double kxxValue = 2.0*Math.PI*ising.R*kxx/ising.L;
					int eta_k_int = Lp*((ky+Lp)%Lp) + (kxx+Lp)%Lp;
					sum += (g_k [(kxx-kx+Lp)%Lp] * ising.findVkSquare(kxxValue,kyValue) + f_k[(kxx-kx+Lp)%Lp])*etaLT_k[eta_k_int];				
					
				}
				change[i] = sum/Lp;
			}
		}
		return change;
	}
	
	double [] simulateLinearWithModDynamics(double [] etaLinear){
		double [] linearTheoryGrowth = new double [Lp*Lp];
		double [] eta_bar = new double[Lp*Lp];
		ising.convolveWithRange(etaLinear, eta_bar, ising.R);
		for (int i = 0; i < Lp*Lp; i++)
			linearTheoryGrowth[i] = -eta_bar[i]*g[i%Lp] + etaLinear[i]*f[i%Lp];
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

	public void findRightSideNoAssumptions(double [] eta_x, double [] config){
		double [] eta_k = new double [Lp*Lp];
		eta_k = fft.calculate2DFT(eta_x);
		coeff.findCoefficientsFromConfig(config);
		for(int i = 0; i < Lp*Lp; i++){
			int currentKx = i%Lp;
			int currentKy = i/Lp;
			double kRxValue = kRChunk*currentKx;
			double kRyValue = kRChunk*currentKy;
			double V_k = ising.findVkSquare(kRxValue, kRyValue);
			double sum = -V_k*eta_k[i];
			for (int j = 0; j < Lp*Lp; j ++){
				int shiftedKx = j%Lp;
				int shiftedKy = j/Lp;
				double coefficient = coeff.aa[shiftedKx][shiftedKy];
				int etaXindex = (currentKx+shiftedKx)%Lp;
				int etaYindex = (currentKy+shiftedKy)%Lp;
				int	etaIndex = etaYindex*Lp + etaXindex;
				sum -= ising.T*coefficient*eta_k[etaIndex];
			}
			rhs[i] = sum;
		}
	}

	public void findRightSide(double [] eta_x, double [] slice){
		double [] eta_k = new double [Lp*Lp];
		eta_k = fft.calculate2DFT(eta_x);
		//coeff.findCoefficientsFromSlice(slice);
		coeff.findCoefficientsFromFile();
		for(int i = 0; i < Lp*Lp; i++){
			//coeff.findCoefficientsFromSlice(slice);
			//coeff.findCoefficientsFromFile();
			int kx = i%Lp;
			int ky = i/Lp;
			double kRxValue = kRChunk*kx;
			double kRyValue = kRChunk*ky;
			double V_k = ising.findVkSquare(kRxValue, kRyValue);
			double sum = -V_k*eta_k[i];
			for (int j = 0; j + kx < Lp; j++)
				sum -= ising.T*coeff.a[j]*eta_k[ i + j ];
			for (int j = Lp - kx; j < Lp; j++)
				sum -= ising.T*coeff.a[j]*eta_k[ i + j - Lp ];
			rhs[i] = sum;
		}
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
		rhs = new double[Lp*Lp];
		eta = new double[Lp*Lp];
		etaLT = new double[Lp*Lp];
		etaLT_k = new double[Lp*Lp];
		phi0 = new double [Lp];
		phi0_bar = new double [Lp];
		g_k = new double[Lp];
		f_k = new double[Lp];
		f = new double[Lp];
		g = new double[Lp];
		findPhi0andPhi0_bar();
	}	

}

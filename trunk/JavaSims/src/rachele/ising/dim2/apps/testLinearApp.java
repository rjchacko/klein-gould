package rachele.ising.dim2.apps;

import java.io.DataInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.awt.Color;
import rachele.ising.dim2.FindCoefficients;
import rachele.ising.dim2.FourierTransformer;
import rachele.ising.dim2.IsingField2D;
import rachele.util.FileUtil;
import scikit.graphics.dim2.Grid;
import scikit.graphics.dim2.Plot;
import scikit.jobs.Control;
import scikit.jobs.Job;
import scikit.jobs.Simulation;
import scikit.jobs.params.ChoiceValue;
import scikit.jobs.params.DoubleValue;
//import scikit.numerics.fft.managed.ComplexDouble2DFFT;
import scikit.dataset.Accumulator;
import scikit.dataset.PointSet;
import scikit.util.*;

public class testLinearApp extends Simulation{
    IsingField2D ising;
    public FindCoefficients coeff;
    public FourierTransformer fft;
    public double [] rhs, etaLT, eta; //right hand side
    public double kRChunk; //=2piR/L
    public Accumulator etaAcc = new Accumulator(1.0);
    public Accumulator etaLTAcc = new Accumulator(1.0);
    public Accumulator etaAcc2 = new Accumulator(1.0);
    public Accumulator etaLTAcc2 = new Accumulator(1.0);
    public Accumulator etaAcc3 = new Accumulator(1.0);
    public Accumulator etaLTAcc3 = new Accumulator(1.0);
    public Accumulator etaAcc4 = new Accumulator(1.0);
    public Accumulator etaLTAcc4 = new Accumulator(1.0);
    public int Lp;
    Grid etaDot = new Grid("ising delta phi");
    Grid etaDotLinear = new Grid("linear theory");
    Grid phiGrid = new Grid("phi field");
    Grid etaDotSF = new Grid("sf phi field");
    Grid etaDotLinearSF = new Grid("sf linear theory");
    Plot slice1 = new Plot("delphi mid slice");
    Plot slice2 = new Plot("rhs mid slice");
    Plot slice3 = new Plot("delphi slice");
    Plot slice4 = new Plot("rhs slice");
    Plot etaVt1 = new Plot("eta v t");
    Plot etaVt2 = new Plot("eta v t LT");
    
	public static void main(String[] args) {
		new Control(new testLinearApp(), "Ising Linear Test");
	}
	public void load(Control c) {
	c.frameTogether("Grids", etaDot, etaDotSF, phiGrid, etaDotLinear, etaDotLinearSF);
	c.frameTogether("Slices", slice1, slice2,slice3, slice4);
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
	}

	public void animate() {
		params.set("Time", ising.time());
		params.set("Mean Phi", ising.mean(ising.phi));
		params.set("Lp", ising.Lp);
		etaDot.registerData(ising.Lp, ising.Lp, ising.delPhi);
		etaDotLinear.registerData(ising.Lp, ising.Lp, rhs);
		phiGrid.registerData(ising.Lp, ising.Lp, ising.phi);
		
		double [] delPhiFT = fft.calculateSF2D(ising.delPhi, false, false);
		double [] rhsFT = fft.calculateSF2D(rhs, false, false);

		//System.out.println(delPhiFT [2*Lp] + " " + delPhiFT [0]);
		double factor = delPhiFT [2*Lp]/DoubleArray.max(delPhiFT);
		double factorRHS = rhsFT[2*Lp]/DoubleArray.max(rhs);
		for (int i = 0; i < Lp; i++){
			delPhiFT[i] *= factor;
			rhsFT[i] *= factorRHS;
		}
		//System.out.println(delPhiFT[0] + " " + delPhiFT[2*Lp] + " " + factor);
		fft.center(delPhiFT);
		fft.center(rhsFT);
		etaDotSF.registerData(ising.Lp, ising.Lp, delPhiFT);
		etaDotLinearSF.registerData(ising.Lp, ising.Lp, rhsFT);
		
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
		slice1.setAutoScale(true);
		slice2.setAutoScale(true);
		slice3.setAutoScale(true);
		slice4.setAutoScale(true);
		
		slice1.registerLines("mid slice", new PointSet(0,1,delPhiMidSlice), Color.BLACK);
		slice2.registerLines("mid slice", new PointSet(0,1,rhsMidSlice), Color.BLUE);
		slice3.registerLines("slice", new PointSet(0,1,delPhiSlice), Color.BLACK);
		slice4.registerLines("slice", new PointSet(0,1,rhsSlice), Color.BLUE);
		
		double [] etaSF = new double[Lp*Lp];
		double [] etaLTSF = new double[Lp*Lp];
		etaSF = fft.calculateSF2D(eta, false, false);
		etaLTSF = fft.calculateSF2D(etaLT, false, false);
		etaLTAcc.accum(ising.time(), etaLTSF[2*Lp]);
		etaAcc.accum(ising.time(), etaSF[2*Lp]);
		etaLTAcc2.accum(ising.time(), etaLTSF[2*Lp + 1]);
		etaAcc2.accum(ising.time(), etaSF[2*Lp + 1]);
		etaLTAcc3.accum(ising.time(), etaLTSF[2*Lp + 2]);
		etaAcc3.accum(ising.time(), etaSF[2*Lp + 2]);
		etaLTAcc4.accum(ising.time(), etaLTSF[2*Lp + 3]);
		etaAcc4.accum(ising.time(), etaSF[2*Lp + 3]);
		etaVt1.registerLines("acc", etaAcc, Color.BLACK);
		etaVt2.registerLines("acc", etaLTAcc, Color.BLACK);
		etaVt1.registerLines("acc2", etaAcc2, Color.BLUE);
		etaVt2.registerLines("acc2", etaLTAcc2, Color.BLUE);
		etaVt1.registerLines("acc3", etaAcc3, Color.RED);
		etaVt2.registerLines("acc3", etaLTAcc3, Color.RED);
		etaVt1.registerLines("acc4", etaAcc4, Color.ORANGE);
		etaVt2.registerLines("acc4", etaLTAcc4, Color.ORANGE);
		
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
		if(params.sget("Init Conditions") == "Read From File")
			readInputParams("../../../research/javaData/configs/inputParams");
		ising = new IsingField2D(params);
		this.Lp = ising.Lp;
		coeff = new FindCoefficients(Lp);
		fft = new FourierTransformer(Lp);
		rhs = new double[Lp*Lp];
		eta = new double[Lp*Lp];
		etaLT = new double[Lp*Lp];
		double [] phiSlice = new double [Lp];
		String fileName = "../../../research/javaData/configs1d/config";
		phiSlice = FileUtil.readConfigFromFile(fileName, Lp);
		double [] phi0 = new double [Lp*Lp];
		for (int i = 0; i < Lp*Lp; i++){
			phi0[i] = phiSlice[i%Lp];
			etaLT[i] = ising.phi[i] - phi0[i];
		}
        while (true) {
    		rhs = ising.simulateLinearWithModDynamics(phi0, etaLT);
    		ising.simulate();
    		for (int i = 0; i < Lp*Lp; i++){
    			etaLT[i] += rhs[i];
    			eta[i] = ising.phi[i] - phi0[i];
    		}
    		Job.animate();
		}	
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
	

}

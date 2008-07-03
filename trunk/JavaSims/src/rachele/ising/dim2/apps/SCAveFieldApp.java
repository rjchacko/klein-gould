package rachele.ising.dim2.apps;

import static scikit.util.Utilities.asList;
import java.awt.Color;
import java.io.DataInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import rachele.ising.dim2.IsingField2D;
import rachele.ising.dim2.StripeClumpFieldSim;
import rachele.util.FileUtil;
import rachele.util.FourierTransformer;
import scikit.dataset.Accumulator;
import scikit.dataset.PointSet;
import scikit.graphics.dim2.Geom2D;
import scikit.graphics.dim2.Grid;
import scikit.graphics.dim2.Plot;
import scikit.jobs.Control;
import scikit.jobs.Job;
import scikit.jobs.Simulation;
import scikit.jobs.params.ChoiceValue;
import scikit.jobs.params.DoubleValue;

/**
* Compares the Ising simulations for square interaction to the linear theory
* following a quench in external field.
* 
* This works for the square shaped interaction shape.  
* Use TestLinearOptApp.java with flexible length scales 
* for circle shape interaction.
* 
* The linear theory is:
* \eta(x,y) = -\int dx' V(x-x',y-y') \eta(x',y')
* 				-T \eta(x,y)/(1-phi_0(x))
* In Fourier space, it becomes a matrix equation.  
* For a given value of ky, we have:
* \etaVec[i] = M[i,j]\etaVec[j]
* where \etaVec is a vector of length Lp; each dimension
* represents an allowed value of kx.
*/
public class SCAveFieldApp extends Simulation{
    IsingField2D ising;
    StripeClumpFieldSim sc;
    FourierTransformer fft;
    //double [] rhs, rhs2D, etaLT, eta, etaLT_k_slice, etaK; //right hand side
    double [] eta, etaK;
    int ky;
    double kRChunk; //=2piR/L
    boolean clearFile;
    
    //RUN OPTIONS
    boolean accEtaValues = true;
    boolean writeToFile = true;
    
    public Accumulator etaAcc;
    public Accumulator etaAcc2;
    public Accumulator etaAcc3;
    public Accumulator etaAcc4;
    public int Lp;
    Grid phiGrid = new Grid("phi field");
    Grid etaDotSF = new Grid("sf phi field");
    Plot etakPlot = new Plot("eta(k)");
    Plot etaVsTimeSim = new Plot("eta v t");
    Plot etaVsTimeLinearK = new Plot("eta v t LT k space");
    Plot etaVsTimeLinear = new Plot("eta v t LT");
    Plot etaVsTimeLC = new Plot("eta v t LC");
    Plot fkPlot = new Plot ("f(k)");
    Plot etakSimPlot = new Plot("eta k sim");
    Plot hSlice = new Plot("Horizontal Slice");
	Plot vSlice = new Plot("Vertical Slice"); 

    
	public static void main(String[] args) {
		new Control(new SCAveFieldApp(), "Ising Linear Test");
	}
	
	public void load(Control c) {
		if(accEtaValues) c.frameTogether("accs", etaVsTimeSim, etaVsTimeLinear, etaVsTimeLinearK, etaVsTimeLC);
		c.frameTogether("Grids", phiGrid, etaDotSF);
		c.frameTogether("Slices", vSlice, hSlice);
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
		//params.addm("dt", 0.01);
		params.addm("J", -1.0);
		params.addm("R", 2000000.0);
		params.addm("Random seed", 0);
		params.add("L/R", 3.0);
		params.add("R/dx", 43.0);
		params.add("kR bin-width", 0.1);
		params.add("Magnetization", 0.0);
		params.addm("ky", 2);
		params.addm("dt", 0.0001);
		params.add("dt new");
		params.add("Time");
		params.add("Mean Phi");
		params.add("Lp");
		flags.add("Clear");
		flags.add("Write 1D Config");
	}

	public void animate() {
		params.set("Time", ising.time());
		params.set("Mean Phi", ising.mean(ising.phi));
		params.set("Lp", ising.Lp);
		phiGrid.registerData(Lp, Lp, ising.phi);
		etaDotSF.registerData(Lp, Lp, etaK);
		etaVsTimeSim.setAutoScale(true);
		etaVsTimeLinear.setAutoScale(true);
		etaVsTimeLinearK.setAutoScale(true);
		etaVsTimeLC.setAutoScale(true);
		hSlice.setAutoScale(true);
		vSlice.setAutoScale(true);
		
		if(accEtaValues){
			etaVsTimeSim.registerLines("acc", etaAcc, Color.BLACK);
			etaVsTimeSim.registerLines("acc2", etaAcc2, Color.BLUE);
			etaVsTimeSim.registerLines("acc3", etaAcc3, Color.RED);
			etaVsTimeSim.registerLines("acc4", etaAcc4, Color.ORANGE);
		}

		double horizontalSlice = params.fget("Horizontal Slice");
		double verticalSlice = params.fget("Vertical Slice");
		
		phiGrid.setDrawables(asList(
				Geom2D.line(0, horizontalSlice, 1, horizontalSlice, Color.GREEN),
				Geom2D.line(verticalSlice, 0, verticalSlice, 1, Color.BLUE)));
		
		hSlice.registerLines("Slice", ising.getHslice(horizontalSlice), Color.GREEN);
		//String fileName = "../../../research/javaData/configs1d/config";
		//double [] phi0 = FileUtil.readConfigFromFile(fileName, ising.Lp);
		hSlice.registerLines("phi0", new PointSet(0, 1, sc.phi0) , Color.BLACK);
		vSlice.registerLines("Slice", ising.getVslice(verticalSlice), Color.BLUE);
		
//		double [] eigenVect = new double [Lp];
//		int i = Lp-1;
//		for (int  j = 0; j < Lp; j++){
//			eigenVect[i+Lp*j]=VV[j][i];
//		}
		
		
		double [] eta_k_slice = new double [Lp];
		double findMax=0.0;
		double findMin=0.0;
		double findMaxe = 0.0;
		double findMine = 0.0;
		for (int i = 0; i < Lp; i++){
			eta_k_slice[i]=etaK[Lp*ky+i];
			if (eta_k_slice[i] > findMax) findMax = eta_k_slice[i];
			if (eta_k_slice[i] < findMin) findMin = eta_k_slice[i];
			if (sc.VV[Lp-1][i] < findMine) findMine = sc.VV[Lp-1][i];
			if (sc.VV[Lp-1][i] > findMaxe) findMaxe = sc.VV[Lp-1][i];
		}
		if(flags.contains("Clear")){
			etaAcc.clear();
			etaAcc2.clear();
			etaAcc3.clear();
			etaAcc4.clear();

			flags.clear();			
		}
		
	}

	public void clear() {
	}
 
	public void run() {
		if (flags.contains("Write 1D Config")){
			write1Dconfig();
			flags.clear();
		}
		clearFile = true;
		ky = params.iget("ky");
		initialize();
		double recordStep = .0001;	
		
//		for (int i = 0; i < Lp*Lp; i++)
//			etaLT[i] = ising.phi[i] - sc.phi0[i%Lp];
//		double [] etaLT_k = new double [Lp*Lp];
//		etaLT_k = fft.calculate2DFT(etaLT);
//		for (int i = ky*Lp; i < Lp*(ky+1); i++)
//			etaLT_k_slice[i-ky*Lp] = etaLT_k[i]; 
//		for (int i = 0; i < Lp*Lp; i++)
//			eta[i] = ising.phi[i] - sc.phi0[i%Lp];

		System.out.println("ky = " + ky);
		calcHspinodal();
		Job.animate();

		while (true) {

			ising.readParams(params);
			//ising.simulateUnstable();
			ising.simulateSimple();
			params.set("dt new", ising.dt);
//			if(accEtaValues){	
//				rhs2D = sc.simulateLinear(etaLT);	
//				//rhs = simulateLinearKbar();
//				rhs = sc.simulateLinearK();
//				sc.etaLinearCombo();
//			}
			for (int i = 0; i < Lp*Lp; i++)
				eta[i] = ising.phi[i] - sc.phi0[i%Lp];
			//etaK = fft.calculate2DFT(eta);			
			etaK = fft.calculate2DSF(eta, false, true);
			if(accEtaValues){

   				etaAcc.accum(ising.time(), Math.pow(etaK[ky*Lp],2));
				etaAcc2.accum(ising.time(),Math.pow(etaK[ky*Lp+1],2));
				etaAcc3.accum(ising.time(), Math.pow(etaK[ky*Lp+2],2));
				etaAcc4.accum(ising.time(), Math.pow(etaK[ky*Lp+3],2)); 				
			}
   		
    		Job.animate();
    		if(writeToFile){
    			if(ising.time() >= recordStep){
    				recordSfDataToFile(etaK);
    				recordStep += .0001;
    			}
    		}
		}	
	}
	
	private void calcHspinodal(){
		double rho = Math.sqrt(1+ising.T/(ising.findVkSquare(IsingField2D.KRsquare,0.0)));
		double hs = rho + (ising.T/2.0)*(Math.log(1.0+rho) - Math.log (1-rho));
		System.out.println("H spinodal for T = " + ising.T + " is " + hs);
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
		//System.out.println("data written for time = " + ising.time());
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
		sc = new StripeClumpFieldSim(ising, ky);
		
		etaAcc = new Accumulator();
		etaAcc2 = new Accumulator();
		etaAcc3 = new Accumulator();
		etaAcc4 = new Accumulator();
		this.Lp = ising.Lp;
		fft = new FourierTransformer(Lp);
		//rhs = new double[Lp];
		//rhs2D = new double[Lp*Lp];
		eta = new double[Lp*Lp];
		etaK = new double [Lp*Lp];
		//etaLT = new double[Lp*Lp];
		//etaLT_k = new double[Lp*Lp];
		//etaLT_k_slice = new double [Lp];
//		phi0 = new double [Lp];
//		phi0_bar = new double [Lp];
		if(params.sget("Init Conditions")=="Read 1D Soln") ising.set1DConfig(sc.phi0);
	}	

	private void write1Dconfig(){
		String configFileName = "../../../research/javaData/configs1d/config";
		FileUtil.deleteFile(configFileName);
		FileUtil.writeConfigToFile(configFileName, ising.Lp, ising.phi);
	}
}

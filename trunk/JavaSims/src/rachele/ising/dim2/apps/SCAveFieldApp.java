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
import scikit.jobs.params.DirectoryValue;
import scikit.jobs.params.DoubleValue;
import scikit.jobs.params.FileValue;

/**
* 
* Compares the Ising simulations for square interaction to the linear theory
* following a quench in external field (or not following a quench depending 
* on initial conditions).
* 
* Plots and averages Structure factor vs time for some max time after the quench.
* Use with svtFieldApp to find good parameters.
* 
* This works for the square shaped interaction shape.  
* Use TestLinearOptApp.java with flexible length scales 
* for circle shape interaction. (Probably not working right now.)
* 
* Seems to be a problem: doesn't work for first run.  Have to reset and restart.  
* Do not know why this is.
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
    
    int accNo = 6;
    Accumulator [] etaAcc = new Accumulator [accNo];
    Accumulator [] etaAccAve = new Accumulator [accNo];
    int [] sfLabel = new int [accNo];
    public int Lp;
    Grid phiGrid = new Grid("phi field");
    Grid etaDotSF = new Grid("sf phi field");
    Plot etakPlot = new Plot("eta(k)");
    Plot etaVsTimeSim = new Plot("eta v t");
    Plot etakSimPlot = new Plot("eta k sim");
    Plot hSlice = new Plot("Horizontal Slice");
	Plot vSlice = new Plot("Vertical Slice"); 

    
	public static void main(String[] args) {
		new Control(new SCAveFieldApp(), "SC Ising Field");
	}
	
	public void load(Control c) {
		c.frameTogether("Everything", phiGrid, etaDotSF, vSlice, hSlice, etaVsTimeSim);
		params.add("1D Input File", new FileValue("/home/erdomi/data/lraim/configs1d/L64R23T0-04h0-8"));
		params.add("Data Dir", new DirectoryValue("/home/erdomi/data/lraim/stripeToClumpInvestigation/ftResults/SCAveFieldApp"));
		params.addm("Zoom", new ChoiceValue("Yes", "No"));
		params.addm("Interaction", new ChoiceValue("Square", "Circle"));
		params.addm("Dynamics?", new ChoiceValue("Langevin No M Convervation"));
		params.add("Init Conditions", new ChoiceValue("Read 1D Soln", "Read From File","Random Gaussian"));
		params.addm("Approx", new ChoiceValue("Slow", "HalfStep", "TimeAdjust", "Phi4","Phi4HalfStep"));
		params.addm("Noise", new DoubleValue(1.0, 0.0, 1.0).withSlider());
		params.addm("Horizontal Slice", new DoubleValue(0.5, 0, 0.9999).withSlider());
		params.addm("Vertical Slice", new DoubleValue(0.5, 0, 0.9999).withSlider());
		params.addm("kR", new DoubleValue(5.135622302, 0.0, 6.0).withSlider());
		params.addm("T", 0.04);
		params.addm("H", 0.80);
		params.addm("tolerance", 0.0001);
		params.addm("J", -1.0);
		params.addm("R", 2000000.0);
		params.addm("Random seed", 0);
		params.add("L/R", 2.7826087);
		params.add("R/dx", 40.0);
		params.add("kR bin-width", 0.1);
		params.add("Magnetization", 0.0);
		params.addm("ky", 2);
		params.addm("Max Time", 100.0);
		params.addm("dt", 0.005);
		params.add("dt new");
		params.add("Time");
		params.add("Mean Phi");
		params.add("Lp");
		params.add("Reps");
	}

	public void animate() {
		params.set("Time", ising.time());
		params.set("Mean Phi", ising.mean(ising.phi));
		params.set("Lp", ising.Lp);
		phiGrid.registerData(Lp, Lp, ising.phi);
		etaDotSF.registerData(Lp, Lp, etaK);
		etaVsTimeSim.setAutoScale(true); etaVsTimeSim.setLogScale(false, true);
		
		
		hSlice.setAutoScale(true);
		vSlice.setAutoScale(true);
		
		for (int i = 0; i < accNo; i ++){
			float colorChunk = (float)i/(float)accNo;
			Color col = Color.getHSBColor(colorChunk, 1.0f, 1.0f);
			StringBuffer sb = new StringBuffer();sb.append("s(t) Ave "); sb.append(i);
			etaVsTimeSim.registerLines(sb.toString(), etaAccAve[i], col);
			StringBuffer sb2 = new StringBuffer();sb2.append("s(t) "); sb2.append(i);
			etaVsTimeSim.registerLines(sb2.toString(), etaAcc[i], col);
		}

		double horizontalSlice = params.fget("Horizontal Slice");
		double verticalSlice = params.fget("Vertical Slice");
		
		phiGrid.setDrawables(asList(
				Geom2D.line(0, horizontalSlice, 1, horizontalSlice, Color.GREEN),
				Geom2D.line(verticalSlice, 0, verticalSlice, 1, Color.BLUE)));
		
		hSlice.registerLines("Slice", ising.getHslice(horizontalSlice), Color.GREEN);
		hSlice.registerLines("phi0", new PointSet(0, 1, sc.phi0) , Color.BLACK);
		vSlice.registerLines("Slice", ising.getVslice(verticalSlice), Color.BLUE);
		
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
			for ( int i = 0; i < accNo; i++)
				etaAcc[i].clear();
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
		int repNo = 0;
		double maxtime = params.fget("Max Time");
		
		while (true) {
			double recordStep = .0001;	
			if(params.sget("Init Conditions")=="Read 1D Soln") ising.set1DConfig(sc.phi0);	
			ising.restartClock();
			for (int i = 0; i < accNo; i ++)
				etaAcc[i].clear();
			//double maxTime = 5.0;
			while (ising.time() < maxtime){
				ising.readParams(params);
				//ising.simulateUnstable();
				ising.simulateSimple();
				//ising.simulate();
				//params.set("dt new", ising.dt);
				if(ising.time() >= recordStep){
					for (int i = 0; i < Lp*Lp; i++)
						eta[i] = ising.phi[i] - sc.phi0[i%Lp];
					etaK = fft.find2DSF(eta, ising.L);
					for (int i = 0; i < accNo; i++){
						etaAcc[i].accum(ising.time(), etaK[sfLabel[i]]);
						//etaAcc[i].accum(ising.time(), etaK[sfLabel[i]+Lp-2*i]);
						etaAccAve[i].accum(ising.time(), etaK[sfLabel[i]]);
						//etaAccAve[i].accum(ising.time(), etaK[sfLabel[i]+Lp-2*i]);
					}

	    			recordStep += 1.0;
	    		}
				//System.out.println("dx = " + ising.dx);
				Job.animate();
			}
			repNo += 1;
			params.set("Reps", repNo);
			recordSfDataToFile();
		}	
	}
	
	private void calcHspinodal(){
		double rho = Math.sqrt(1+ising.T/(ising.findVkSquare(IsingField2D.KRsquare,0.0)));
		double hs = rho + (ising.T/2.0)*(Math.log(1.0+rho) - Math.log (1-rho));
		System.out.println("H spinodal for T = " + ising.T + " is " + hs);
	}
	
	public void recordSfDataToFile(){
		
		String message1 = "#Field theory average of SF vs time for several values of k.  Stripe to clump. ";
		String fileName = params.sget("Data Dir") + File.separator + "f0";
		StringBuffer fileBuffer = new StringBuffer(); fileBuffer.append(fileName);
		for (int i=0; i < accNo; i ++){
			System.out.println("start " + i);
			StringBuffer mb = new StringBuffer();
			mb.append("# sf label = ");	mb.append(sfLabel[i]); mb.append(" kR value = ");
			double krvalue = 2*ising.R*Math.PI*(sfLabel[i])/ising.L;
			mb.append(krvalue);			
			fileBuffer.deleteCharAt(fileBuffer.length()-1);	fileBuffer.append(i); fileName = fileBuffer.toString();
			String message2 = mb.toString();
			FileUtil.initFile(fileName, params, message1, message2);		
			FileUtil.printAccumToFile(fileName, etaAccAve[i]);
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
		sc = new StripeClumpFieldSim(ising, ky, params.sget("1D Input File"),params.sget("Data Dir"));
		
		for (int i = 0; i < accNo; i++){
			etaAcc[i] = new Accumulator();
			etaAccAve[i] = new Accumulator (); etaAccAve[i].enableErrorBars(true);
			sfLabel[i] = ky*Lp+i;
		}	
		this.Lp = ising.Lp;
		fft = new FourierTransformer(Lp);
		eta = new double[Lp*Lp];
		etaK = new double [Lp*Lp];
		if(params.sget("Init Conditions")=="Read 1D Soln") ising.set1DConfig(sc.phi0);
	}	

	
	private void write1Dconfig(){
		String configFileName = "../../../research/javaData/configs1d/config";
		FileUtil.deleteFile(configFileName);
		FileUtil.writeConfigToFile(configFileName, ising.Lp, ising.phi);
	}
}

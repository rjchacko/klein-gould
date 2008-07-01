
package rachele.ising.dim2.apps;


import static scikit.util.Utilities.*;
import scikit.dataset.Accumulator;
import scikit.dataset.PointSet;
import java.awt.Color;
import java.io.*;
import rachele.ising.dim2.ConjugateGradientMin;
import rachele.ising.dim2.IsingField2D;
import rachele.ising.dim2.SteepestDescentMin;
import rachele.ising.dim2.StructureFactor;
import scikit.graphics.dim2.Geom2D;
import scikit.graphics.dim2.Grid;
import scikit.graphics.dim2.Plot;
import scikit.jobs.Control;
import scikit.jobs.Job;
import scikit.jobs.Simulation;
import scikit.jobs.params.ChoiceValue;
import scikit.jobs.params.DoubleValue;
import rachele.util.FileUtil;
import static java.lang.Math.pow;


public class IsingField2DApp extends Simulation {
    Grid grid = new Grid("Phi(x)");
    Grid grid2 = new Grid("Phi(x)");
    Grid sfGrid = new Grid("S(k)");
    Grid delPhiGrid = new Grid("DelPhi");
	Plot hSlice = new Plot("Horizontal Slice");
	Plot vSlice = new Plot("Vertical Slice");    
	Plot slicePlot1 = new Plot("ky=0 Slice");
	Plot slicePlot2 = new Plot("ky=k0 Slice");
	Plot sfHor = new Plot("H SF");
	Plot sfVert = new Plot("V SF");
	Plot structurePeakV = new Plot("Ver Structure Factor");
	Plot structurePeakH = new Plot("Hor Structure factor");
	Plot sfPeakBoth = new Plot("Both Structure factors");
	Plot freeEnergyPlot = new Plot("Free Energy");
	Plot freeEnergyTempPlot = new Plot("Free Energy vs Temp");
	Plot landscape = new Plot("Free Energy Landscape");
	Plot brLandscape = new Plot("Br Free Energy Landscape");
	Plot ring = new Plot("Circle SF Ring");
	Plot ringInput = new Plot("Input for Ring");
	Plot dS_dtPlot = new Plot("SF change");
	StructureFactor sf;
    IsingField2D ising;
    SteepestDescentMin opt;
    ConjugateGradientMin min;
    boolean cgInitialized = false;
	boolean initFile = false;
    Accumulator landscapeFiller;
    Accumulator brLandscapeFiller;
    Accumulator sfChange;
    Accumulator sfSliceAcc;
    public double [] inputSlice;
    public int lastClear;
    public int maxi=0;
    
	public static void main(String[] args) {
		new Control(new IsingField2DApp(), "Ising Field");
	}
	
	public void load(Control c) {
		//uncomment next line
		c.frameTogether("Grids", grid, delPhiGrid, sfGrid, freeEnergyPlot);
		c.frameTogether("Slices", vSlice, hSlice, slicePlot1, slicePlot2);
		c.frame(hSlice);
		//structurePeakH, freeEnergyPlot, sfPeakBoth, sfHor, sfVert);
		//c.frameTogether("SF", sfHor, sfVert);
		params.addm("Zoom", new ChoiceValue("Yes", "No"));
		params.addm("Interaction", new ChoiceValue("Square", "Circle"));
		params.addm("Dynamics?", new ChoiceValue("Langevin No M Convervation", "Langevin Conserve M","Conjugate Gradient Min", 
				"Steepest Decent"));
		params.add("Init Conditions", new ChoiceValue("Read From File","Random Gaussian", 
				 "Artificial Stripe 3", "Artificial Stripe 2","Constant", "Read 1D Soln"));
		params.addm("Approx", new ChoiceValue("Slow", "HalfStep", "TimeAdjust", "Phi4","Phi4HalfStep"));
		//params.addm("Plot FEvT", new ChoiceValue("Off", "On"));
		params.addm("Noise", new DoubleValue(0.0, 0.0, 1.0).withSlider());
		params.addm("Stripe Strength", new DoubleValue(0.0, 0.0, 1.0).withSlider());
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
		params.add("Free Energy");
		flags.add("Write Config");
		flags.add("Write 1D Config");
		flags.add("Clear");
		//flags.add("SF");
		//flags.add("new RS");
		//flags.add("Integrate");
		landscapeFiller = new Accumulator(.01);
		brLandscapeFiller = new Accumulator(.01);
		sfChange = new Accumulator(1);
		sfSliceAcc = new Accumulator(.1);
	}
	
	public void animate() {

		ising.readParams(params);
		
		if (params.sget("Zoom").equals("Yes")) {
			grid.setAutoScale();
			grid2.setAutoScale();
			delPhiGrid.setAutoScale();
		}
		else {
			grid.setScale(-1, 1);
			grid2.setScale(-1, 1);
			delPhiGrid.setScale(0, 1);
		}
		
		freeEnergyPlot.setAutoScale(true);
		slicePlot1.setAutoScale(true);
		slicePlot2.setAutoScale(true);
		hSlice.setAutoScale(true);
		vSlice.setAutoScale(true);
		structurePeakV.setAutoScale(true);
		structurePeakH.setAutoScale(true);
		sfPeakBoth.setAutoScale(true);
		landscape.setAutoScale(true);
		brLandscape.setAutoScale(true);
		sfVert.setAutoScale(true);
		sfHor.setAutoScale(true);
		ring.setAutoScale(true);
		dS_dtPlot.setAutoScale(true);
		
		params.set("Free Energy", ising.freeEnergy);
//		for(int i=0; i < ising.Lp; i++){
//			int j = ising.Lp*ising.Lp/2 + i;
//			sf.sFactor[j]=0;
//		}
		sfGrid.registerData(ising.Lp, ising.Lp, sf.sFactor);
		grid.registerData(ising.Lp, ising.Lp, ising.phi);
		grid2.registerData(ising.Lp, ising.Lp, ising.phi);
		//String dyn = params.sget("Dynamics?");
		landscape.registerLines("FE landscape", landscapeFiller, Color.BLACK);
		brLandscape.registerLines("FE landscape", brLandscapeFiller, Color.BLUE);
		delPhiGrid.registerData(ising.Lp, ising.Lp, ising.phiVector);
		
		
		slicePlot1.registerLines("Slice", sf.get_sfSlice(), Color.RED);
		slicePlot2.registerLines("Slice", getInput(), Color.YELLOW);

		for (int i =0; i<ising.Lp; i++){
			int j = (ising.Lp*ising.Lp)/2 +2*ising.Lp + i;
			sfChange.accum(i, sf.sFactor[j]);
		}
		dS_dtPlot.registerLines("sf change slice", sfChange, Color.BLACK);

		double horizontalSlice = params.fget("Horizontal Slice");
		double verticalSlice = params.fget("Vertical Slice");
		
		grid.setDrawables(asList(
				Geom2D.line(0, horizontalSlice, 1, horizontalSlice, Color.GREEN),
				Geom2D.line(verticalSlice, 0, verticalSlice, 1, Color.BLUE)));

		delPhiGrid.setDrawables(asList(
				Geom2D.line(0, horizontalSlice, 1, horizontalSlice, Color.RED),
				Geom2D.line(verticalSlice, 0, verticalSlice, 1, Color.BLACK)));
		
		hSlice.registerLines("Slice", ising.getHslice(horizontalSlice), Color.GREEN);
		//String fileName = "../../../research/javaData/configs1d/config";
		//double [] phi0 = FileUtil.readConfigFromFile(fileName, ising.Lp);
		//hSlice.registerLines("phi0", new PointSet(0, 1, phi0) , Color.BLACK);
		vSlice.registerLines("Slice", ising.getVslice(verticalSlice), Color.BLUE);
		freeEnergyPlot.registerLines("Free Energy", ising.getFreeEnergyAcc(), Color.MAGENTA);
		
		if(ising.circleInt() == true){
			ring.registerLines("RING", sf.getRingFT(), Color.black);
			//ringInput.registerLines("Input", sf.getRingInput(), Color.black);
			structurePeakV.registerLines("Peak Value", sf.getPeakC(), Color.BLACK);
		}else{
			structurePeakV.registerLines("Vertical Peak", sf.getPeakV(), Color.CYAN);
			sfPeakBoth.registerLines("Hortizontal Peak", sf.getPeakH(), Color.ORANGE);
			sfPeakBoth.registerLines("Vertical Peak", sf.getPeakV(), Color.CYAN);
			structurePeakH.registerLines("Horizontal Peak", sf.getPeakH(), Color.ORANGE);
			sfHor.registerLines("Hor SF Ave", sf.getAccumulatorHA(), Color.BLACK);
			sfHor.registerLines("Hor SF", sf.getAccumulatorH(), Color.ORANGE);
			sfVert.registerLines("Vert SF Ave", sf.getAccumulatorVA(), Color.BLACK);
			sfVert.registerLines("Vert SF", sf.getAccumulatorV(), Color.CYAN);
		}	
 
		freeEnergyTempPlot.registerLines("Stripe FE", ising.getStripeFreeEnergyAcc(), Color.GREEN);
		freeEnergyTempPlot.registerLines("Clump FE", ising.getClumpFreeEnergyAcc(),Color.PINK);
		freeEnergyTempPlot.registerLines("fe", ising.getEitherFreeEnergyAcc(), Color.BLACK);
		
		if (flags.contains("Clear")){// || lastClear > 1000) {
			ising.getFreeEnergyAcc().clear();
			sf.getPeakH().clear();
			sf.getPeakV().clear();
			sf.getPeakC().clear();
			sf.getPeakHslope().clear();
			sf.getPeakVslope().clear();
			sf.getAccumulatorVA().clear();
			sf.getAccumulatorHA().clear();
			ising.aveCount = 0;
			lastClear = 0;
		}
		if(flags.contains("SF")){
			for (int i = 0; i < 5000; i ++){
				ising.simulate();
			}
			maxi=sf.clumpsOrStripes(ising.phi);
		}
		if(flags.contains("new RS")) ising.randomizeSeed(params.iget("Random seed", 0));
		if (flags.contains("Integrate")){
			double V_k = -0.2067483214;
			double a0H = integrate(true);
			double a0V = integrate(false);
			double slope1 =2*(-V_k-ising.T*a0H);
			double slope2 =2*(-V_k-ising.T*a0V);
			double linearSlope = 2*(-V_k-ising.T/(1.0-pow(ising.mean(ising.phi),2)));
			System.out.println("slope1 = " + slope1 + " slope2 = " + slope2 + " linear slope = " + linearSlope);
		}
		if(flags.contains("Write 1D Config")) write1Dconfig();
		flags.clear();
	}
	
	public void clear() {
		cgInitialized = false;
		initFile = false;
	}
	
	public PointSet getInput(){
		return new PointSet(0, 1, inputSlice);
	}
	
	public void run() {
		if(params.sget("Init Conditions") == "Read From File")
			readInputParams("../../../research/javaData/configs/inputParams");
		ising = new IsingField2D(params);
		inputSlice = new double [ising.Lp];
		double binWidth = params.fget("kR bin-width");
		//binWidth = ising.KRcircle / floor(IsingField2D.KR_SP/binWidth);
        sf = new StructureFactor(ising.Lp, ising.L, ising.R, binWidth, ising.dt);
		sf.setBounds(0.1, 14);
		int recordSteps = 0;
		maxi=sf.clumpsOrStripes(ising.phi);
		
        while (true) {
        	ising.readParams(params);
        	if (flags.contains("Write Config"))	writeConfiguration();
			params.set("Time", ising.time());
			params.set("Mean Phi", ising.mean(ising.phi));
			if(params.sget("Dynamics?") == "Conjugate Gradient Min"){
				if(cgInitialized == false){
					min.initialize();   
					System.out.println("CG initialized");
					cgInitialized = true;					
				}
				min.step(ising.t);
				ising.accFreeEnergy.accum(ising.t, min.freeEnergy);
				landscapeFiller = min.getLandscape();					
				brLandscapeFiller = min.getBracketLandscape();
				ising.t += 1;
			}else if(params.sget("Dynamics?") == "Steepest Decent"){
				opt.step();
				ising.t += 1;
				ising.accFreeEnergy.accum(ising.t, opt.freeEnergy);
				cgInitialized = false;
				landscapeFiller = opt.getLandscape();
				brLandscapeFiller = opt.getBracketLandscape();
			}else{
				cgInitialized = false;
				ising.simulate();
				//ising.simulateUnstable();
			}
			sf.getAccumulatorV().clear();
			sf.getAccumulatorH().clear();
			double [] input = new double [ising.Lp*ising.Lp];
			for (int i = 0; i < ising.Lp*ising.Lp; i ++)
				input[i] = 1.0/(1.0-Math.pow(ising.phi[i], 2));
			for (int i = 0; i < ising.Lp; i ++)
				inputSlice[i] = input[i];
			sf.accumulateAll(ising.time(), input);
			
			//sf.accumulateAll(ising.time(), ising.coarseGrained());
			
			if (ising.time() > recordSteps){
				//sf.accumMin(ising.coarseGrained(), params.fget("kR"));
				boolean circleOn=false;
				sf.accumulateMelt(circleOn, ising.phi, maxi);
				//sf.accumulateAll(ising.t, ising.delPhi);

				sf.accumulateAll(ising.t, ising.phi);
				recordSFvTime();
				//record3Ddata();
				recordSteps += 0.1;
				writeDataToFile();
			}
			Job.animate();
//			int half = ising.Lp*ising.Lp/2 + ising.Lp/2;
//			double ratio = sf.sFactor[half + 2]/sf.sFactor[half + 4];
//			System.out.println("ratio = " + ratio);
			sfChange.clear();
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
	
	public void writeConfiguration(){
		String configFileName = "../../../research/javaData/configs/inputConfig";
		String inputFileName = "../../../research/javaData/configs/inputParams";
		FileUtil.deleteFile(configFileName);
		FileUtil.deleteFile(inputFileName);
		writeInputParams(inputFileName);	
		writeConfigToFile(configFileName);
		System.out.println("Config writtern to file");
	}
	
	public void writeConfigToFile(String FileName){
		try {
			File pathFile = new File(FileName);
			DataOutputStream dos = new DataOutputStream(new FileOutputStream(pathFile, true));
			for (int i = 0; i < ising.Lp*ising.Lp; i ++){
				dos.writeInt(i);
				dos.writeChar('\t');
				dos.writeDouble(ising.phi[i]);
				dos.writeChar('\n');
			}
			dos.close();
		} catch (IOException ex){
			ex.printStackTrace();
		}
	}
	
	public void writeInputParams(String FileName){
		try {
			File inputFile = new File(FileName);
			DataOutputStream dos = new DataOutputStream(new FileOutputStream(inputFile, true));
			
			dos.writeDouble(params.fget("J"));
			dos.writeChar('\t');
//			dos.writeDouble(params.fget("H"));
//			dos.writeChar('\t');
			dos.writeDouble(params.fget("R"));
			dos.writeChar('\t');
			dos.writeDouble(params.fget("L/R"));
			dos.writeChar('\t');
			dos.writeDouble(params.fget("R/dx"));
			dos.writeChar('\t');
			dos.writeChar('\n');
			dos.close();
		}catch(IOException ex){
			ex.printStackTrace();
		}
	}

	public void initFile(String file, boolean SvH){
		FileUtil.deleteFile(file);
		if(SvH){
			FileUtil.printlnToFile(file, " # SF vs H data ");			
			FileUtil.printlnToFile(file, " # Temperature = ", ising.T);
			FileUtil.printlnToFile(file, " # Data = H, S(k*), Free Energy, time");
		}else{
			FileUtil.printlnToFile(file, " # SF vs T data ");				
			//FileUtil.printlnToFile(file, " # External field = ", ising.H);
			FileUtil.printlnToFile(file, " # Data = H, S(k*), Free Energy, time");
		}
		FileUtil.printlnToFile(file, " # Density = ", ising.DENSITY);		
	}
	
	public void initSFvTimeFile(String file){
		FileUtil.deleteFile(file);
		FileUtil.printlnToFile(file, " # SF vs time data ");			
		FileUtil.printlnToFile(file, " # Temperature = ", ising.T);
		//FileUtil.printlnToFile(file, " # H", ising.H);
		FileUtil.printlnToFile(file, " # Density = ", ising.DENSITY);		
	}
	
	public void recordSFvTime(){
		if (params.sget("Interaction")=="Square"){
			String dataFileV = "../../../research/javaData/sfData/sfv";
			String dataFileH = "../../../research/javaData/sfData/sfh";
			if (initFile == false){
				initSFvTimeFile(dataFileV);
				initSFvTimeFile(dataFileH);
				initFile = true;
			}
			FileUtil.printlnToFile(dataFileH, ising.time(), sf.peakValueH());
			FileUtil.printlnToFile(dataFileV, ising.time(), sf.peakValueV());					
			System.out.println("Data written to file for time = " + ising.time());
		}else{
			System.out.println("no write to file for non-square yet");
		}
	}
	
	public void writeDataToFile(){
		boolean SvH = true;
		if (params.sget("Interaction")=="Square"){
				String dataFileV = "../../../research/javaData/sfData/sV";
				String dataFileH = "../../../research/javaData/sfData/sH";
				if (initFile == false){
					initFile(dataFileV, SvH);
					initFile(dataFileH, SvH);
					initFile = true;
				}
				if (SvH){
					FileUtil.printlnToFile(dataFileH, params.fget("H"), sf.peakValueH(), ising.freeEnergy, ising.time());
					FileUtil.printlnToFile(dataFileV, params.fget("H"), sf.peakValueV(), ising.freeEnergy, ising.time());					
				}else{
					FileUtil.printlnToFile(dataFileH, ising.T, sf.peakValueH(), ising.freeEnergy, ising.time());
					FileUtil.printlnToFile(dataFileV, ising.T, sf.peakValueV(), ising.freeEnergy, ising.time());
				}			
				System.out.println("Data written to file for time = " + ising.time());
		}else if(params.sget("Interaction")== "Circle"){
			String dataStripe = "../../../research/javaData/sfData/dataStripe";
			String dataClump = "../../../research/javaData/sfData/dataClump";
			if (initFile == false){
				initFile(dataStripe, SvH);
				initFile(dataClump, SvH);
				initFile = true;
			}
			if(SvH){
				FileUtil.printlnToFile(dataClump, params.fget("H"), sf.peakValueC(), ising.freeEnergy, ising.time());
				FileUtil.printlnToFile(dataStripe, params.fget("H"), sf.peakValueS(), ising.freeEnergy, ising.time());					
			}else{
				FileUtil.printlnToFile(dataClump, ising.T, sf.peakValueC(), ising.freeEnergy, ising.time());
				FileUtil.printlnToFile(dataStripe, ising.T, sf.peakValueS(), ising.freeEnergy, ising.time());
			}
			System.out.println("Data written to file for time = " + ising.time());
		}
	}

//	private void record3Ddata() {
//			String dataFile = "../../../research/javaData/sfData/3d";
//
//			if (initFile == false){
//				initFile(dataFile, true);
//				initFile = true;
//			}
//			for(int i = 0; i < ising.Lp; i ++){
//				FileUtil.printlnToFile(dataFile, ising.time(), (double)i, sf.sFactor[ising.Lp*ising.Lp/2 + i]);
//			}
//	}
	
	private double integrate(boolean direction){
		double integralSum = 0;
		if (direction){
			for(int k = 0; k < ising.Lp; k ++){
				double lineSum = 0;
				for(int i = 0; i < ising.Lp; i ++){
					int point = ising.Lp*k + i;
					lineSum += 1.0/(1-pow(ising.phi[point],2));
				}
				double sumAve = lineSum/(double)ising.Lp;
				integralSum += sumAve;
			}
		}else{
			for (int k = 0; k < ising.Lp; k ++){
				double lineSum = 0;
				for(int i = 0; i < ising.Lp; i ++){
					int point = i*ising.Lp + k;
					lineSum += 1.0/(1-pow(ising.phi[point],2));
				}
				double sumAve = lineSum/(double)ising.Lp;
				integralSum += sumAve;
			}

		}
		double ave = integralSum/(double)ising.Lp;
		return ave;
	}
	
	private void write1Dconfig(){
		String configFileName = "../../../research/javaData/configs1d/config";
		FileUtil.deleteFile(configFileName);
		double[] slice = new double [ising.Lp];
		for (int i = 0; i < ising.Lp; i ++)
			slice[i] = ising.phi[i];
		FileUtil.writeConfigToFile(configFileName, ising.Lp, slice);
	}

}


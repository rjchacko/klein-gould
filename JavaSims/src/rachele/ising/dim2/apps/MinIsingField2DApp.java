package rachele.ising.dim2.apps;

import static scikit.util.Utilities.asList;
import java.awt.Color;
import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import rachele.ising.dim2.IsingField2D;
import rachele.util.FileUtil;
import rachele.util.FourierTransformer;
import scikit.dataset.PointSet;
import scikit.graphics.dim2.Geom2D;
import scikit.graphics.dim2.Grid;
import scikit.graphics.dim2.Plot;
import scikit.jobs.Control;
import scikit.jobs.Job;
import scikit.jobs.Simulation;
import scikit.jobs.params.ChoiceValue;
import scikit.jobs.params.DoubleValue;

public class MinIsingField2DApp extends Simulation{
    Grid grid = new Grid("Phi(x)");
    Grid sfGrid = new Grid("S(k)");
    Grid delPhiGrid = new Grid("DelPhi");
	Plot hSlice = new Plot("Horizontal Slice");
	Plot vSlice = new Plot("Vertical Slice");    
	
	FourierTransformer fft;
    IsingField2D ising;
   	boolean initFile = false;
    public int lastClear;
    public int maxi=0;
    double [] sf;
    
	public static void main(String[] args) {
		new Control(new MinIsingField2DApp(), "Ising Field");
	}
	
	public void load(Control c) {
		c.frameTogether("Grids", grid, sfGrid, vSlice, hSlice);
		params.addm("Zoom", new ChoiceValue("Yes", "No"));
		params.addm("Interaction", new ChoiceValue("Square", "Circle"));
		params.addm("Dynamics?", new ChoiceValue("Langevin No M Convervation", "Langevin Conserve M","Conjugate Gradient Min", 
				"Steepest Decent"));
		params.add("Init Conditions", new ChoiceValue("Read From File","Random Gaussian", 
				 "Artificial Stripe 3", "Artificial Stripe 2","Constant", "Read 1D Soln"));
		params.addm("Approx", new ChoiceValue("Slow", "HalfStep", "TimeAdjust", "Phi4","Phi4HalfStep"));
		params.addm("Noise", new DoubleValue(0.0, 0.0, 1.0).withSlider());
		params.addm("Stripe Strength", new DoubleValue(0.0, 0.0, 1.0).withSlider());
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
		params.add("Free Energy");
		flags.add("Write Config");
		flags.add("Write 1D Config");
		flags.add("Clear");

	}
	
	public void animate() {

		ising.readParams(params);
		if (params.sget("Zoom").equals("Yes")) {
			grid.setAutoScale();
			delPhiGrid.setAutoScale();
		}
		else {
			grid.setScale(-1, 1);
			delPhiGrid.setScale(0, 1);
		}
		
		hSlice.setAutoScale(true);
		vSlice.setAutoScale(true);
		params.set("Free Energy", ising.freeEnergy);
		sfGrid.registerData(ising.Lp, ising.Lp, sf);
		grid.registerData(ising.Lp, ising.Lp, ising.phi);
		delPhiGrid.registerData(ising.Lp, ising.Lp, ising.phiVector);
		double horizontalSlice = params.fget("Horizontal Slice");
		double verticalSlice = params.fget("Vertical Slice");
		
		grid.setDrawables(asList(
				Geom2D.line(0, horizontalSlice, 1, horizontalSlice, Color.GREEN),
				Geom2D.line(verticalSlice, 0, verticalSlice, 1, Color.BLUE)));

		delPhiGrid.setDrawables(asList(
				Geom2D.line(0, horizontalSlice, 1, horizontalSlice, Color.RED),
				Geom2D.line(verticalSlice, 0, verticalSlice, 1, Color.BLACK)));
		
		hSlice.registerLines("Slice", ising.getHslice(horizontalSlice), Color.GREEN);
		String fileName = "../../../research/javaData/configs1d/config";
		double [] phi0 = FileUtil.readConfigFromFile(fileName, ising.Lp);
		hSlice.registerLines("phi0", new PointSet(0, 1, phi0) , Color.BLACK);
		vSlice.registerLines("Slice", ising.getVslice(verticalSlice), Color.BLUE);

		if (flags.contains("Clear")){// || lastClear > 1000) {
			ising.getFreeEnergyAcc().clear();
			ising.aveCount = 0;
			lastClear = 0;
		}
		if(flags.contains("new RS")) ising.randomizeSeed(params.iget("Random seed", 0));
		if(flags.contains("Write 1D Config")) write1Dconfig();
		flags.clear();
	}
	
	public void clear() {
		initFile = false;
	}
	
	public void run() {
		if(params.sget("Init Conditions") == "Read From File")
			readInputParams("../../../research/javaData/configs/inputParams");
		ising = new IsingField2D(params);
		sf = new double [ising.Lp*ising.Lp];
		fft = new FourierTransformer(ising.Lp);
		int recordSteps = 0;

        while (true) {
        	ising.readParams(params);
        	if (flags.contains("Write Config"))	writeConfiguration();
			params.set("Time", ising.time());
			params.set("Mean Phi", ising.mean(ising.phi));
			ising.simulate();
			//sf = fft.calculate2DSF(ising.phi, true, true);
			if (ising.time() > recordSteps){
				sf = fft.calculate2DSF(ising.phi, false, false);
				recordSteps += 1;
				recordSFvTime();
			}
			Job.animate();
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
	
	private void writeConfiguration(){
		String configFileName = "../../../research/javaData/configs/inputConfig";
		String inputFileName = "../../../research/javaData/configs/inputParams";
		FileUtil.deleteFile(configFileName);
		FileUtil.deleteFile(inputFileName);
		writeInputParams(inputFileName);	
		writeConfigToFile(configFileName);
		System.out.println("Config writtern to file");
	}
	
	private void writeConfigToFile(String FileName){
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
	
	private void writeInputParams(String FileName){
		try {
			File inputFile = new File(FileName);
			DataOutputStream dos = new DataOutputStream(new FileOutputStream(inputFile, true));
			
			dos.writeDouble(params.fget("J"));
			dos.writeChar('\t');
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


	
	public void initSFvTimeFile(String file){
		FileUtil.deleteFile(file);
		FileUtil.printlnToFile(file, " # SF vs time data ");			
		FileUtil.printlnToFile(file, " # Temperature = ", ising.T);
		FileUtil.printlnToFile(file, " # Density = ", ising.DENSITY);		
	}
	
	public void recordSFvTime(){
		String dataFile0 = "../../../research/javaData/sfData/s0";
		String dataFile1v = "../../../research/javaData/sfData/s1v";
		String dataFile1h = "../../../research/javaData/sfData/s1h";
		String dataFile2v = "../../../research/javaData/sfData/s2v";
		String dataFile2h = "../../../research/javaData/sfData/s2h";
		if (initFile == false){
			initSFvTimeFile(dataFile0);
			initSFvTimeFile(dataFile1v);
			initSFvTimeFile(dataFile1h);
			initSFvTimeFile(dataFile2v);
			initSFvTimeFile(dataFile2h);
			initFile = true;
		}
		FileUtil.printlnToFile(dataFile0, ising.time(), sf[0]);
		FileUtil.printlnToFile(dataFile1h, ising.time(), sf[1]);
		FileUtil.printlnToFile(dataFile2h, ising.time(), sf[2]);
		FileUtil.printlnToFile(dataFile1v, ising.time(), sf[ising.Lp]);
		FileUtil.printlnToFile(dataFile2v, ising.time(), sf[ising.Lp*2]);

		System.out.println("Data written to file for time = " + ising.time());
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


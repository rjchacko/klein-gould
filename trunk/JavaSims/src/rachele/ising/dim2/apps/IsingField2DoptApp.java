package rachele.ising.dim2.apps;


import static java.lang.Math.floor;
//import static scikit.util.Utilities.*;
import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.EOFException;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import rachele.ising.dim2.IsingField2Dopt;
import rachele.util.FileUtil;
import rachele.ising.dim2.StructureFactorOpt;
import scikit.graphics.dim2.Grid;
import scikit.jobs.Control;
import scikit.jobs.Job;
import scikit.jobs.Simulation;
import scikit.dataset.Accumulator;
import scikit.jobs.params.ChoiceValue;
import scikit.jobs.params.DoubleValue;
import scikit.graphics.dim2.Plot;
import java.awt.Color;

public class IsingField2DoptApp extends Simulation{
    Grid grid = new Grid("Phi(x)");
    Grid delPhiGrid = new Grid("delPhi(x)");
    Grid sfGrid = new Grid("S(k)");
    Plot fePlot = new Plot("Free Energy");
	StructureFactorOpt sf;
	IsingField2Dopt ising;
	boolean initFile = false;
	boolean showFE = true;
	Accumulator freeEnergy;
    
    public static void main(String[] args) {
		new Control(new IsingField2DoptApp(), "Ising Field");
	}
	
	public void load(Control c){
		c.frameTogether("Grids", grid, sfGrid, delPhiGrid);
		if (showFE) c.frame(fePlot);
		params.addm("Zoom", new ChoiceValue("Yes", "No"));
		params.addm("Interaction", new ChoiceValue( "Circle","Square" ));
		params.addm("Theory", new ChoiceValue("Exact", "Slow Near Edge", "Dynamic dt"));
		params.addm("Dynamics?", new ChoiceValue("Langevin No M Convervation", "Langevin Conserve M"));
		params.add("Init Conditions", new ChoiceValue("Random Gaussian", "Read From File"));
		params.addm("Noise", new DoubleValue(0, 0, 1.0).withSlider());
		params.addm("T", 0.02);
		params.addm("H", 0.8);
		params.addm("Rx", 2490000.0);
		params.addm("Ry", 2160000.0);
		params.add("L", 6000000.0);
		params.add("dx", 60000.0);
		params.add("Random seed", 0);
		params.add("Magnetization", 0.0);
		params.addm("range change", 0.01);
		params.add("dt", 1.0);
		params.add("Time");
		params.add("Free Energy");
		params.add("Pot");
		params.add("Ent");
		flags.add("Write Config");
		flags.add("Record FE");
		flags.add("Clear");
		flags.add("Write Config");

	}
	
	public void animate() {
		ising.readParams(params);
		if (params.sget("Zoom").equals("Yes"))grid.setAutoScale();
		else grid.setScale(-1, 1);
		fePlot.setAutoScale(true);
		if(showFE){
			freeEnergy.accum(ising.t, ising.freeEnergy);
			fePlot.registerLines("FE", freeEnergy, Color.RED);
		}
		sfGrid.registerData(ising.Lp, ising.Lp, sf.sFactor);
		grid.registerData(ising.Lp, ising.Lp, ising.phi);
		delPhiGrid.registerData(ising.Lp, ising.Lp, ising.phiVector);
		params.set("Rx", ising.Rx);
		params.set("Ry", ising.Ry);
		params.set("Free Energy", ising.freeEnergy);
		params.set("Pot", ising.potAccum);
		params.set("Ent", ising.entAccum);
		params.set("dt", ising.dt);
		if(flags.contains("Clear")) freeEnergy.clear();
		if(flags.contains("Record FE")) recordTvsFE();
		if(ising.recordTvsFE){
			recordTvsFE();
			double newT = ising.T - 0.001;
			params.set("T", newT);			
			ising.recordTvsFE = false;
		}else if(ising.recordHvsFE){
			recordHvsFE();
			double newH = ising.H + 0.005;
			params.set("H", newH);			
			ising.recordHvsFE = false;			
		}
		flags.clear();
	}
	
	public void clear() {
		initFile = false;
	}
	
	public void run() {
		boolean recordSFtoFile = true;
		ising = new IsingField2Dopt(params);
		if(params.sget("Init Conditions") == "Read From File") readInitialConfiguration();
		double binWidth = 0.1;
		binWidth = IsingField2Dopt.KR_SP / floor(IsingField2Dopt.KR_SP/binWidth);
        sf = new StructureFactorOpt(ising.Lp, ising.L);
        if(showFE) freeEnergy = new Accumulator(1.0);
		int steps = 1;
		int recordSteps = 10;
//		if (ising.t < 500.0){
//			params.set("Time", ising.time());
//			ising.simulate();
//			ising.adjustRanges();			
//		}
//		steps = 500;
        while (true) {
        	ising.readParams(params);
        	if (flags.contains("Write Config"))	writeConfiguration();
        	params.set("Time", ising.time());
			ising.simulate();
			//ising.adjustRanges();
			if (ising.t > steps){
				sf.takeFT(ising.phi);
				Job.animate();
				steps += 1;
			}
			if(recordSFtoFile){	
				if (ising.time() > recordSteps){
					sf.shiftSFactor(); //This is necessary to do once before peak value finds
					double peakValueH = sf.findSquarePeak(ising.Rx, true);
					double peakValueV = sf.findSquarePeak(ising.Ry, false);
					recordSFvTime(peakValueH, peakValueV);
					recordSteps += 10;
				}
			}
        }
 	}
	
	public void recordSFvTime(double peakValueH, double peakValueV){
		if (params.sget("Interaction")=="Square"){
			String dataFileV = "../../../research/javaData/sfData/sfv";
			String dataFileH = "../../../research/javaData/sfData/sfh";
			if (initFile == false){
				initSFvTimeFile(dataFileV);
				initSFvTimeFile(dataFileH);
				initFile = true;
			}
			FileUtil.printlnToFile(dataFileH, ising.time(), peakValueH);
			FileUtil.printlnToFile(dataFileV, ising.time(), peakValueV);					
			System.out.println("Data written to file for time = " + ising.time());
		}else{
			System.out.println("no write to file for non-square yet");
		}
	}
	
	public void initSFvTimeFile(String file){
		FileUtil.deleteFile(file);
		FileUtil.printlnToFile(file, " # SF vs time data ");			
		FileUtil.printlnToFile(file, " # Temperature = ", ising.T);
		FileUtil.printlnToFile(file, " # H", ising.H);
		FileUtil.printlnToFile(file, " # Density = ", ising.DENSITY);		
	}
	
	public void recordTvsFE(){
		String file = "../../../research/javaData/feData/fe";
		FileUtil.printlnToFile(file, ising.T, ising.freeEnergy);
		System.out.println("Wrote to file: T = " + ising.T + " FE = " + ising.freeEnergy);
	}

	public void recordHvsFE(){
		String file = "../../../research/javaData/feData/fe";
		FileUtil.printlnToFile(file, ising.H, ising.freeEnergy);
		System.out.println("Wrote to file: H = " + ising.H + " FE = " + ising.freeEnergy);
	}
	

	public void initFile(String file, boolean SvH){
		FileUtil.deleteFile(file);
		if(SvH){
			FileUtil.printlnToFile(file, " # SF vs H data ");			
			FileUtil.printlnToFile(file, " # Temperature = ", ising.T);
			FileUtil.printlnToFile(file, " # Data = H, S(k*), Free Energy, time");
		}else{
			FileUtil.printlnToFile(file, " # SF vs T data ");				
			FileUtil.printlnToFile(file, " # External field = ", ising.H);
			FileUtil.printlnToFile(file, " # Data = H, S(k*), Free Energy, time");
		}
		FileUtil.printlnToFile(file, " # Density = ", ising.DENSITY);		
	}
	public void writeConfiguration(){
		String configFileName = "../../../research/javaData/configs/inputConfigOpt";
		String inputFileName = "../../../research/javaData/configs/inputParamsOpt";
		FileUtil.deleteFile(configFileName);
		FileUtil.deleteFile(inputFileName);
		writeInputParams(inputFileName);	
		writeConfigToFile(configFileName);
		System.out.println("Config writtern to file");
	}
	
	public void writeInputParams(String FileName){
		try {
			File inputFile = new File(FileName);
			DataOutputStream dos = new DataOutputStream(new FileOutputStream(inputFile, true));
			
			dos.writeDouble(params.fget("H"));
			dos.writeChar('\t');
			dos.writeDouble(params.fget("Rx"));
			dos.writeChar('\t');
			dos.writeDouble(params.fget("Ry"));
			dos.writeChar('\t');
			dos.writeDouble(params.fget("L"));
			dos.writeChar('\t');
			dos.writeDouble(params.fget("dx"));
			dos.writeChar('\t');
			dos.writeChar('\n');
			dos.close();
		}catch(IOException ex){
			ex.printStackTrace();
		}
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
	public void readInitialConfiguration(){
		try{
			File myFile = new File("../../../research/javaData/configs/inputConfigOpt");
			DataInputStream dis = new DataInputStream(new FileInputStream(myFile));
			int spaceIndex;
			double phiValue;
			try{
				while(true){
					spaceIndex =dis.readInt();
					dis.readChar();       // throws out the tab
					phiValue = dis.readDouble();
					dis.readChar();
					ising.phi[spaceIndex] = phiValue;
				}
			} catch (EOFException e) {
			}

		} catch (FileNotFoundException e) {
			System.err.println("FileStreamsTest: " + e);
		} catch (Exception ex) {
			ex.printStackTrace();
		}
	}
}

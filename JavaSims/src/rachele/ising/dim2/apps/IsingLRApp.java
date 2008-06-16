package rachele.ising.dim2.apps;

import java.io.DataInputStream;
import java.io.EOFException;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;

import rachele.ising.dim2.*;
import scikit.graphics.dim2.Grid;
import scikit.jobs.Control;
import scikit.jobs.Job;
import scikit.jobs.Simulation;
import scikit.jobs.params.ChoiceValue;
import static scikit.util.Utilities.format;
import rachele.util.*;


public class IsingLRApp extends Simulation {
	
	Grid grid = new Grid("Long Range Ising Model");
	Grid sfGrid = new Grid("SF");
	int dx;
	IsingLR sim;
	public FourierTransformer fft;
	double [] sFactor;
    boolean clearFile;
    boolean takeAverages = true;
   	double [] sfTimeArray;
	
	public static void main(String[] args) {
		new Control(new IsingLRApp(), "Ising Model");
	}
	
	public void load(Control c) {
		c.frameTogether("grids", grid, sfGrid);		
		params.addm("Dynamics", new ChoiceValue("Ising Glauber","Kawasaki Glauber", "Kawasaki Metropolis",  "Ising Metropolis"));
		params.addm("init", new ChoiceValue("Read From File", "Random"));
		params.add("Random seed", 0);
		params.add("L", 1<<8);
		params.add("R", 86);
		params.add("Initial magnetization", 0.0);
		params.addm("T", 0.04);
		params.addm("J", -1.0);
		params.addm("h", 0.8);
		params.addm("dt", 0.1);
		params.add("time");
		params.add("magnetization");
		params.add("Lp");
		params.add("Reps");
		flags.add("Write Config");
	}	
	
	public void animate() {
		params.set("time", format(sim.time()));
		params.set("magnetization", format(sim.magnetization()));
		sim.setParameters(params);
		params.set("Lp", sim.L/dx);
		
		grid.registerData(sim.L/dx, sim.L/dx, sim.getField(dx));
		sfGrid.registerData(sim.L/dx, sim.L/dx, sFactor);
		if(flags.contains("Write Config")) writeConfigToFile();
		flags.clear();
	}
	
	public void clear() {
	}
	
	public void run() {
		initialize();
		fft = new FourierTransformer((int)(sim.L/dx));
		sFactor = new double [sim.L/dx*sim.L/dx];
		double step = 0.10;
		if(takeAverages){
			double maxTime = 10.0;
			int ky = 2;
			int sfLabel = ky*sim.L/dx + 0;
			sfTimeArray = new double[(int)((maxTime/step)+1)];
			int repNo = 0;
			while (true) {
				readInitialConfiguration();
				sim.restartClock();
				int recordInt = 0;
				int recordStep = 0;
				while (sim.time() < maxTime){
					sim.step();
					Job.animate();
					if (sim.time() > recordStep){
						sFactor = fft.calculate2DSF(sim.getField(dx), false, false);
						sfTimeArray[recordInt] += sFactor[sfLabel];
						recordStep += step;
						recordInt +=1;
					}
				}	
				repNo += 1;
				params.set("Reps", repNo);
				writeArray(repNo, step);
			}
		}else{
			int recordStep = 0;
			while (true) {
				sim.step();
				if (sim.time() > recordStep){
					sFactor = fft.calculate2DSF(sim.getField(dx), false, false);
					recordSfDataToFile(sFactor);
					recordStep += step;
				}
				Job.animate();
			}
		}
	}

	public void writeArray(int repNo, double stepSize){

		String fileName = "../../../research/javaData/stripeToClumpInvestigation/monteCarloData/sf";
		FileUtil.deleteFile(fileName);
		int length = sfTimeArray.length;
		for (int i = 0; i < length; i++){
			double value = sfTimeArray[i]/repNo;
			double time = stepSize*(double)i;
			FileUtil.printlnToFile(fileName, time, value);
		}
	}
	
	public void initialize(){
		sim = new IsingLR(params);
		sim.setField(params.fget("Initial magnetization"));
		dx = 1;
		if(params.sget("init") == "Read From File") readInitialConfiguration();
		clearFile = true;
	}
	
	public void writeConfigToFile(){
		String configFileName = "../../../research/javaData/stripeToClumpInvestigation/monteCarloData/configs/config";
		FileUtil.deleteFile(configFileName);
		FileUtil.writeConfigToFile(configFileName, (sim.L/dx)*(sim.L/dx), sim.getField(dx));
	}
	
	public void readInitialConfiguration(){
		try{
			File myFile = new File("../../../research/javaData/stripeToClumpInvestigation/monteCarloData/configs/config");
			DataInputStream dis = new DataInputStream(new FileInputStream(myFile));
			int spaceIndex;
			double phiValue;
			int Lp = sim.L/dx;
			try{
				while(true){
					spaceIndex =dis.readInt();
					dis.readChar();       // throws out the tab
					phiValue = dis.readDouble();
					dis.readChar();
					sim.spins.set(spaceIndex%Lp,spaceIndex/Lp,(int)phiValue);
//					[spaceIndex] = phiValue;
				}
			} catch (EOFException e) {
			}

		} catch (FileNotFoundException e) {
			System.err.println("FileStreamsTest: " + e);
		} catch (Exception ex) {
			ex.printStackTrace();
		}
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
		int ky = 2;
		int Lp = sim.L/dx;
		FileUtil.printlnToFile(file0, sim.time(), data[ky*Lp]*data[ky*Lp]);
		FileUtil.printlnToFile(file1, sim.time(), data[ky*Lp+1]*data[ky*Lp+1]);
		FileUtil.printlnToFile(file2, sim.time(), data[ky*Lp+2]*data[ky*Lp+2]);
		FileUtil.printlnToFile(file3, sim.time(), data[ky*Lp+3]*data[ky*Lp+3]);
		FileUtil.printlnToFile(file4, sim.time(), data[ky*Lp+4]*data[ky*Lp+4]);
		FileUtil.printlnToFile(file5, sim.time(), data[ky*Lp+5]*data[ky*Lp+5]);
		//System.out.println("data written for time = " + ising.time());
	}	
	
}

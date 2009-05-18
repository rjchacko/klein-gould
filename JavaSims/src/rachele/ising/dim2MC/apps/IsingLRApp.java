package rachele.ising.dim2MC.apps;

import java.awt.Color;
import java.io.DataInputStream;
import java.io.EOFException;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
//import rachele.ising.dim2.*;
import rachele.ising.dim2MC.IsingLR;
import scikit.dataset.Accumulator;
import scikit.graphics.dim2.Grid;
import scikit.graphics.dim2.Plot;
import scikit.jobs.Control;
import scikit.jobs.Job;
import scikit.jobs.Simulation;
import scikit.jobs.params.ChoiceValue;
import static scikit.util.Utilities.format;
import rachele.util.*;
import static java.lang.Math.*;


public class IsingLRApp extends Simulation {
	
	Grid grid = new Grid("Long Range Ising Model");
	Grid sfGrid = new Grid("SF");
	Plot sfPlot = new Plot("sf plot");
	int dx;
	IsingLR sim;
	public FourierTransformer fft;
	double [] sFactor;
	Accumulator sf_k;
	Accumulator sfTheoryAcc, sfTheory2Acc;
	Accumulator sfAcc, sfAccv;
    boolean clearFile;
    boolean takeAverages;
   	double [] sfTimeArray;
	
	public static void main(String[] args) {
		new Control(new IsingLRApp(), "Ising Model");
	}
	
	public void load(Control c) {
		c.frameTogether("grids", grid, sfGrid, sfPlot);		
		params.addm("Dynamics", new ChoiceValue("Ising Glauber","Kawasaki Glauber", "Kawasaki Metropolis",  "Ising Metropolis"));
		params.addm("init", new ChoiceValue( "Random", "Read From File"));
		params.add("Random seed", 0);
		params.add("L", 1<<7);
		params.add("R", 46);//1<<6);
		params.add("Initial magnetization", 0.0);
		params.addm("T", 0.02);
		params.addm("J", -1.0);
		params.addm("h", 0.543);
		params.addm("dt", 1/(double)(1<<3));
		params.add("time");
		params.add("magnetization");
		params.add("Lp");
		params.add("Reps");
		flags.add("Clear");
		sf_k = new Accumulator();
		sfTheoryAcc = new Accumulator();
		sfAcc = new Accumulator();
		sfAccv = new Accumulator();
		
	}	
	
	public void animate() {
		params.set("time", format(sim.time()));
		params.set("magnetization", format(sim.magnetization()));
		sim.setParameters(params);
		params.set("Lp", sim.L/dx);
		
		grid.registerData(sim.L/dx, sim.L/dx, sim.getField(dx));
		sfGrid.registerData(sim.L/dx, sim.L/dx, sFactor);
		sfPlot.registerPoints("sf", sfAcc, Color.RED);
		sfPlot.registerPoints("sfv", sfAccv, Color.BLUE);
		if(flags.contains("Write Config")) writeConfigToFile();
		if(flags.contains("Clear")){
			sfAccv.clear();
			sfAcc.clear();
		}
		flags.clear();
	}
	
	public void clear() {
		sfAcc.clear();
		sfAccv.clear();
	}
	
	public void run() {
		initialize();
		fft = new FourierTransformer((int)(sim.L/dx));
		sFactor = new double [sim.L/dx*sim.L/dx];
		double step = 0.10;
		int recordStep = 0;
		while (true) {
			sim.step();
			if (sim.time() > recordStep){
				sFactor = fft.calculate2DSF(sim.getField(dx), false, false);
//				recordSfDataToFile(sFactor);
				sfAcc.accum(sim.time(), sFactor[2]);
				sfAccv.accum(sim.time(), sFactor[sim.L*2]);
				recordStep += step;
			}
			Job.animate();
		}
	}

	public void structureTheoryAcc(double time, double kRmax){
		System.out.println("time = " + time);
		double kRbinWidth = 0.2;
		sf_k = new Accumulator(kRbinWidth);
		sfTheoryAcc = new Accumulator(kRbinWidth);
		sfTheory2Acc = new Accumulator(kRbinWidth);

		for(int i = 0; i < sim.L/(2*dx); i++){
			double kR=2*Math.PI*(double)i*sim.R/sim.L;
			if (kR < kRmax){
				double pot = (kR == 0) ? 1 : Math.sin(kR)/kR; 
				double D = -pot/sim.T-1;
				//double V = sim.L/dx;
				double sf = (exp(2*time*D)*(1 + 1/D)-1/D);
	//			exp(M*D*t)*(1 + 1/D) - 1/D
				//System.out.println("time = " + time + "sf = " + sf);
				sfTheoryAcc.accum(kR,sf);
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
		//sim.setField(params.fget("Initial magnetization"));
		sim.randomizeField(params.fget("Initial magnetization"));		
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
		String file0 = "../../../research/javaData/stripeToClumpInvestigation/monteCarloData/s0";
		String file1 = "../../../research/javaData/stripeToClumpInvestigation/monteCarloData/s1";
		String file2 = "../../../research/javaData/stripeToClumpInvestigation/monteCarloData/s2";
		String file3 = "../../../research/javaData/stripeToClumpInvestigation/monteCarloData/s3";
		String file4 = "../../../research/javaData/stripeToClumpInvestigation/monteCarloData/s4";
		String file5 = "../../../research/javaData/stripeToClumpInvestigation/monteCarloData/s5";
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

package rachele.ising.dim2.apps;

import java.awt.Color;
import java.io.DataInputStream;
import java.io.EOFException;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import rachele.ising.dim2.*;
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

public class MonteCarloDataApp extends Simulation{

	Grid grid = new Grid("Long Range Ising Model");
	Grid sfGrid = new Grid("SF");
	Plot sfkPlot = new Plot("sf_k plot");
	Plot sftPlot = new Plot("sf_t plot");
	int dx;
	IsingLR sim;
	public FourierTransformer fft;
	double [] sFactor;
	Accumulator sf_kAcc, sf_tAcc;
	Accumulator sf_kTheoryAcc, sf_kTheory2Acc, sf_tTheoryAcc, sf_tTheory2Acc;
    boolean clearFile;
    String averages;
   	double [] sfTimeArray;
	
	public static void main(String[] args) {
		new Control(new MonteCarloDataApp(), "Monte Carlo");
	}
	
	public void load(Control c) {

		//averages = "S_k_DO";//take s(k) averages for disorder to order case
		//averages = "S_t_DO";//take s(t) averages for disorder to order case
		//averages = "S_k_SC;;//take s(k) averages for stripe to clump case
		averages = "S_t_SC";//take s(k) averages for stripe to clump case
		//averages = "None";
		
		c.frameTogether("grids", grid, sfGrid);		
		if (averages == "S_k_DO"){
			c.frame(sfkPlot);
			sf_kAcc = new Accumulator();
		}else if(averages == "S_t_DO"){
			c.frame(sftPlot);
			sf_tAcc = new Accumulator();
		}else if (averages == "S_k_SC"){
			c.frame(sfkPlot);
			sf_kAcc = new Accumulator();
		}else if (averages == "S_t_SC"){
			c.frame(sftPlot);
			sf_tAcc = new Accumulator();
		}
		params.addm("Dynamics", new ChoiceValue("Ising Glauber","Kawasaki Glauber", "Kawasaki Metropolis",  "Ising Metropolis"));
		params.addm("init", new ChoiceValue( "Random", "Read From File"));
		params.add("Random seed", 0);
		params.add("L", 1<<9);
		params.add("R", 1<<6);
		params.add("Initial magnetization", 0.0);
		params.addm("T", 0.04);
		params.addm("J", -1.0);
		params.addm("h", 0.0);
		params.addm("dt", 1/(double)(1<<3));
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
		if (averages == "S_k_DO"){
			sfkPlot.registerLines("sf(k)", sf_kAcc, Color.BLACK);
			sfkPlot.registerLines("theory", sf_kTheoryAcc, Color.BLUE);
		}else if(averages == "S_t_DO"){
			sftPlot.registerLines("st(k)", sf_tAcc, Color.BLACK);
			sftPlot.registerLines("theory", sf_tTheoryAcc, Color.BLUE);			
		}
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
		if(averages == "S_t_DO"){
			double maxTime = 2.0;
			int sfLabel = sf_t_theory(maxTime);
			sfTimeArray = new double[(int)((maxTime/step)+1)];
			int repNo = 0;
			while (true) {
				if(params.sget("init") == "Read From File") readInitialConfiguration();
				else sim.randomizeField(params.fget("Initial magnetization"));
				sim.restartClock();
				int recordInt = 0;
				int recordStep = 0;
				while (sim.time() < maxTime){
					sim.step();
					Job.animate();
					if (sim.time() > recordStep){
						sFactor = fft.calculate2DSF(sim.getField(dx), false, false);
						sfTimeArray[recordInt] += sFactor[sfLabel];
						sf_tAcc.accum(sim.time(),sFactor[sfLabel]);
						recordStep += step;
						recordInt +=1;
					}
				}	
				repNo += 1;
				params.set("Reps", repNo);
				writeArray(repNo, step);
			}
		}else if (averages == "S_k_DO"){
			double maxTime = 0.2;
			double kRmax= 25;
			sf_k_theory(0.25, kRmax);
			System.out.println("take data time = " + maxTime);
			int repNo = 0;
			sfTimeArray = new double[sim.L/dx];			
			while (true) {
				if(params.sget("init") == "Read From File") readInitialConfiguration();
				else sim.randomizeField(params.fget("Initial magnetization"));
				sim.restartClock();
				while (sim.time() < maxTime){
					sim.step();
					Job.animate();
				}	
				sFactor = fft.calculate2DSF(sim.getField(dx), false, false);
				
				for (int i = 1; i < sim.L/(2*dx); i++ ){
					System.out.println("print time = " + sim.time());
					double kRValue = 2*Math.PI*sim.R*(double)i/sim.L;
					if (kRValue < kRmax) sf_kAcc.accum(kRValue, sFactor[i]);
				}
				repNo += 1;
				params.set("Reps", repNo);
			}
		}else if(averages == "S_t_SC"){
			double maxTime = 2.0;//max time after quench time
			int sfLabel = findBestkR();
			int sfLabelVert = sfLabel*sim.L/dx;
			sfTimeArray = new double[(int)((maxTime/step)+1)];
			int repNo = 0;
			double initTime = 15.0;
			while (true) {
				sim.randomizeField(0.0);
				sim.restartClock();
				boolean init = true;

				while(init){
					while(sim.time() < initTime){ 
						sim.step();
						System.out.println("init time = " + sim.time());
					}
					sFactor = fft.calculate2DSF(sim.getField(dx), false, false);
					if(sFactor[sfLabel]>sFactor[sfLabelVert]) init = false;
				}
				sim.restartClock();
				int recordInt = 0;
				int recordStep = 0;
				while (sim.time() < maxTime){
					sim.step();
					Job.animate();
					if (sim.time() > recordStep){
						sFactor = fft.calculate2DSF(sim.getField(dx), false, false);
						sfTimeArray[recordInt] += sFactor[sfLabel];
						sf_tAcc.accum(sim.time(),sFactor[sfLabel]);
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

	public int sf_t_theory(double maxTime){
		sf_tAcc.clear();
		sf_tAcc = new Accumulator(sim.dTime());
		sf_tTheoryAcc = new Accumulator(sim.dTime());
		sf_tTheory2Acc = new Accumulator(sim.dTime());
		int kRint = findBestkR();
		double kRValue = 2*PI*kRint*sim.R/sim.L;
		while (sim.time() < maxTime){
			sim.step();
			//Job.animate();
			System.out.println("theory time = " + sim.time() + " maxtime = " + maxTime + " kr value = " + kRValue);
			double sf = theoryPoint(kRValue, sim.time());
			sf_tTheoryAcc.accum(sim.time(), sf);
		}
		return kRint;
	}
	
	public int findBestkR(){
		int kRint = (int)(IsingLR.kRpeak*sim.L/(2*Math.PI*sim.R));
		double trialkR1 = 2*PI*kRint*sim.R/sim.L;
		double trialkR2 = 2*PI*(kRint+1)*sim.R/sim.L;
		if (abs(IsingLR.kRpeak-trialkR1) > abs(IsingLR.kRpeak-trialkR2)){
			kRint += 1;
		}		
		return kRint;
	}
	
	public void sf_k_theory(double time, double kRmax){
		sf_kAcc.clear();
		System.out.println("time = " + time);
		double kRbinWidth = 0.2;
		sf_kAcc = new Accumulator(kRbinWidth);
		sf_kTheoryAcc = new Accumulator(kRbinWidth);
		sf_kTheory2Acc = new Accumulator(kRbinWidth);

		for(int i = 0; i < sim.L/(2*dx); i++){
			double kR=2*Math.PI*(double)i*sim.R/sim.L;
			if (kR < kRmax){
				double sf = theoryPoint(kR,time);
				sf_kTheoryAcc.accum(kR,sf);
			}
		}
	}
	
	public double theoryPoint(double kR, double time){
		double pot = (kR == 0) ? 1 : Math.sin(kR)/kR; 
		double D = -pot/sim.T-1;
		double sf = (exp(2*time*D)*(1 + 1/D)-1/D);		
		return sf;
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

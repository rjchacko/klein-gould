package rachele.ising.dim2MC.apps;

import static java.lang.Math.exp;
import static scikit.util.Utilities.format;

import java.awt.Color;
import java.io.DataInputStream;
import java.io.EOFException;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;

import rachele.ising.dim2MC.IsingLR;
import rachele.util.FileUtil;
import rachele.util.FourierTransformer;
import scikit.dataset.Accumulator;
import scikit.graphics.dim2.Grid;
import scikit.graphics.dim2.Plot;
import scikit.jobs.Control;
import scikit.jobs.Job;
import scikit.jobs.Simulation;
import scikit.jobs.params.ChoiceValue;

public class IsingKawJumpApp extends Simulation{
	
	Grid grid = new Grid("Long Range Ising Model");
	Grid sfGrid = new Grid("SF");
	Plot sfPlot = new Plot("Sfs");

	IsingLR sim;
	public FourierTransformer fft;
	
	int dx;
	double [] sFactor;
    boolean clearFile;
    boolean takeAverages;
   	double [] sfTimeArray;
   	
	Accumulator sf_k, sf_tv, sf_th, sf_th2, sf_tv2, sf_tv3, sf_th3;
	Accumulator sfTheoryAcc, sfTheory2Acc;
   	
	public static void main(String[] args) {
		new Control(new IsingKawJumpApp(), "Ising Model");
	}
	
	public void load(Control c) {
		c.frameTogether("grids", grid, sfGrid, sfPlot);		
		params.addm("Dynamics", new ChoiceValue("Kawasaki Glauber", "Kawasaki Metropolis", "Ising Glauber",  "Ising Metropolis"));
		params.addm("init", new ChoiceValue( "Random", "Read From File"));
		params.add("Random seed", 0);
		params.add("L", 1<<7);
		params.add("R", 46);
		params.add("Initial magnetization", 0.0);
		params.addm("T", 0.01);
		params.addm("J", -1.0);
		params.addm("h", 0.0);
		params.addm("dt", 1.0);//1/(double)(1<<3));
		params.addm("Jump Range", 1);
		params.addm("kx", 2);
		params.addm("ky", 2);
		params.add("time");
		flags.add("Write Config");
		flags.add("Clear");
		sf_k = new Accumulator();
		sfTheoryAcc = new Accumulator();
		sf_th = new Accumulator();
		sf_tv = new Accumulator();
		sf_th2 = new Accumulator();
		sf_tv2 = new Accumulator();	
		sf_th3 = new Accumulator();
		sf_tv3 = new Accumulator();	
	}	
	
	public void animate() {
		params.set("time", format(sim.time()));
		sim.setParameters(params);
		
		sfPlot.setAutoScale(true);
		grid.registerData(sim.L/dx, sim.L/dx, sim.getField(dx));
		sfGrid.registerData(sim.L/dx, sim.L/dx, sFactor);
		sfPlot.registerPoints("sf_h", sf_th, Color.BLUE);
		sfPlot.registerPoints("sf_v", sf_tv, Color.RED);
		sfPlot.registerPoints("sf_h2", sf_th2, Color.GREEN);
		sfPlot.registerPoints("sf_v2", sf_tv2, Color.ORANGE);
		sfPlot.registerPoints("sf_h3", sf_th3, Color.PINK);
		sfPlot.registerPoints("sf_v3", sf_tv3, Color.CYAN);
		
		if(flags.contains("Write Config")) writeConfigToFile();
		if(flags.contains("Clear")){
			sf_tv.clear();
			sf_th.clear();			
		}
		flags.clear();
	}
	
	public void clear() {
		sf_tv.clear();
		sf_th.clear();
		sf_th2.clear();
		sf_tv2.clear();
		sf_th3.clear();
		sf_tv3.clear();
	}
	
	public void run() {
		initialize();
		fft = new FourierTransformer((int)(sim.L/dx));
		sFactor = new double [sim.L/dx*sim.L/dx];
		double step = 0.10;
		int recordStep = 0;
		int kx = params.iget("kx");
		int ky = params.iget("ky");
		int k2 = 5;
		int k3 = 10;
		while (true) {
			sim.step();
			if (sim.time() > recordStep){
				sFactor = fft.calculate2DSF(sim.getField(dx), false, false);
				sf_th.accum(sim.time(), sFactor[kx]);
				sf_th.accum(sim.time(), sFactor[sim.L-kx]);
				sf_tv.accum(sim.time(), sFactor[ky*sim.L]);
				sf_tv.accum(sim.time(), sFactor[(sim.L-ky)*sim.L]);
				sf_th2.accum(sim.time(), sFactor[k2]);
				sf_th2.accum(sim.time(), sFactor[sim.L-k2]);
				sf_tv2.accum(sim.time(), sFactor[k2*sim.L]);
				sf_tv2.accum(sim.time(), sFactor[(sim.L-k2)*sim.L]);
				sf_th3.accum(sim.time(), sFactor[k3]);
				sf_th3.accum(sim.time(), sFactor[sim.L-k3]);
				sf_tv3.accum(sim.time(), sFactor[k3*sim.L]);
				sf_tv3.accum(sim.time(), sFactor[(sim.L-k3)*sim.L]);
				//				recordSfDataToFile(sFactor);
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
		sim.jumpRange = params.iget("Jump Range");
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

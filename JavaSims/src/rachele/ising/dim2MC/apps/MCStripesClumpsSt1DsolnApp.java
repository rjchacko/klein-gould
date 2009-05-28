package rachele.ising.dim2MC.apps;

import static java.lang.Math.PI;
import static java.lang.Math.abs;
import static scikit.util.Utilities.format;

import java.awt.Color;
import java.io.File;
import java.util.Random;

import rachele.ising.dim2MC.IsingLR;
import rachele.util.FileUtil;
import rachele.util.FourierTransformer;
import scikit.dataset.Accumulator;
import scikit.dataset.PointSet;
import scikit.graphics.dim2.Grid;
import scikit.graphics.dim2.Plot;
import scikit.jobs.Control;
import scikit.jobs.Job;
import scikit.jobs.Simulation;
import scikit.jobs.params.ChoiceValue;
import scikit.jobs.params.DirectoryValue;
import scikit.jobs.params.FileValue;
//import static scikit.util.DoubleArray.*;

/**
* 
* Monte Carlo Simulation to produce structure factor vs time averages
* for stripe to clump transitions.  Averages various structure factor modes
* for a given ky value and several kx values.
* 
*/
public class MCStripesClumpsSt1DsolnApp extends Simulation{

	Grid grid = new Grid("Long Range Ising Model");
	Grid sfGrid = new Grid("Structure Factor");
	Plot sftPlot = new Plot("sf_t plot");
	Plot slicePlot = new Plot("Slice");
	int dx;
	IsingLR sim;
	public FourierTransformer fft;
	double [] sFactor;
	Random random = new Random();
	Accumulator sf_kAcc;
	int accNo = 20;
	double [] inputSlice;
	double [] unstableSlice;
	Accumulator [] sf_tAveAcc = new Accumulator [accNo];
	Accumulator [] sf_tAcc = new Accumulator [accNo];
	int [] sfLabel = new int [accNo];
	int ky;
	
	public static void main(String[] args) {
		new Control(new MCStripesClumpsSt1DsolnApp(), "Stripes to Clumps");
	}

	public void load(Control c) {

		c.frameTogether("Data",grid,sfGrid,sftPlot,slicePlot);
		params.add("Data Dir",new DirectoryValue("/Users/erdomi/data/lraim/stripeToClumpInvestigation/mcResults/conservedOP_SC/testruns"));
		params.add("Input 1D File",new FileValue("/Users/erdomi/data/lraim/configs1dAutoName/L256R91T0.08h0.65"));
		params.add("Input 1D Unstable File",new FileValue("/Users/erdomi/data/lraim/configs1dAutoName/L256R91T0.02h0.543"));
		params.addm("Dynamics", new ChoiceValue("Ising Glauber", "Kawasaki Glauber",  "Kawasaki Metropolis",  "Ising Metropolis"));
		params.add("Random seed", 0);
		params.add("L", 128);//1<<9);
		params.add("R", 46);//1<<6);
		params.add("Jump Range", 1);
		params.add("Initial magnetization", 0.0);
//		params.addm("T", 0.096548);
//		params.addm("T", 0.07952444);
		params.addm("T", 0.02);
		params.addm("J", -1.0);
		params.addm("h", 0.6);
		params.addm("dt", 1.0);//1/(double)(1<<1));
		params.addm("maxTime", 50.0);
		params.addm("ky", 2);
		params.addm("kx int", 1);
		params.add("time");
//		params.add("magnetization");
		params.add("mean phi");
		params.add("Lp");
		params.add("Reps");
		flags.add("Write Config");
	}	

	public void animate() {
		sftPlot.setLogScale(false, true);
		params.set("time", format(sim.time()));
//		params.set("magnetization", format(sim.magnetization()));
		double meanPhi = (double)(sim.spins.sumAll())/(sim.L*sim.L);
		params.set("mean phi", meanPhi);
		sim.setParameters(params);
		params.set("Lp", sim.L/dx);
		grid.registerData(sim.L/dx, sim.L/dx, sim.getField(dx));
//		for(int i = 0; i < sim.L; i++)
//			sFactor[i] = 0;
//		sfGrid.registerData(sim.L/dx, sim.L/dx,sFactor);
		slicePlot.registerLines("Slice", sim.getAveHslice(), Color.BLACK);
		slicePlot.registerLines("Slice Input",  new PointSet(0, 1, inputSlice), Color.GREEN);
		slicePlot.registerLines("Slice Unstable",  new PointSet(0, 1, unstableSlice), Color.RED);
		for (int i = 0; i < accNo; i ++){
			StringBuffer sb = new StringBuffer();sb.append("s(t) Ave "); sb.append(i);
			float colorChunk = (float)i/(float)accNo;
			Color col = Color.getHSBColor(colorChunk, 1.0f, 1.0f);
			sftPlot.registerLines(sb.toString(), sf_tAveAcc[i], col);
			//sb = new StringBuffer();sb.append("s(t) "); sb.append(i);
			//sftPlot.registerLines(sb.toString(), sf_tAcc[i], col);
		}
		//sftPlot.registerLines("Blah", sf_tAcc[0], Color.BLACK);
		
		if(flags.contains("Write Config")) writeConfigToFile();
		flags.clear();
	}

	public void clear() {
	}

	public void run() {
		initialize();
		fft = new FourierTransformer((int)(sim.L/dx));
		sFactor = new double [sim.L/dx*sim.L/dx];
		double step = 0.20;
		double maxTime = params.fget("maxTime");//max time after quench time
		//double quenchH = 0.9;
		int repNo = 0;
//		int sfLabel = findBestkR();
		System.out.println(sfLabel);

		while (true) {
			for (int i = 0; i < accNo; i ++)
				sf_tAcc[i].clear();
			initializeStripes();
//			sim.randomizeField(params.fget("Initial magnetization"));
			sim.restartClock();
			sFactor = fft.calculate2DSF(sim.getField(dx), false, true);
			sim.restartClock();
			int recordInt = 0;
			int recordStep = 0;
			step = 0.25;
			while (sim.time() < maxTime){
				sim.step();
				Job.animate();
				if (sim.time() > recordStep){
//					sFactor = fft.find2DSF(sim.getField(dx), sim.L);
					sFactor = fft.calculate2DSF(sim.getField(dx), false, true);
					collect(sFactor,sfLabel);
					recordStep += step;
					recordInt +=1;
				}
			}	
			repNo += 1;
			params.set("Reps", repNo);
			writeStSCtoFile();
		}
	}

	private void collect(double [] sFactor, int [] sfLabelHor){
		for (int i = 0; i < accNo; i ++){
			sf_tAveAcc[i].accum(sim.time(),sFactor[sfLabel[i]]);			
			sf_tAcc[i].accum(sim.time(),sFactor[sfLabel[i]]);
		}
	}
	
	private double kRvalue(int label){
		double kRvalue = 2*PI*label*sim.R/(sim.L); 
		return kRvalue;
	}

	private void initializeStripes(){
		
		//now load the Ising lattice with the proper probability
		for (int i = 0; i < sim.L*sim.L; i++){
			double prob = (inputSlice[i%sim.L]+1.0)/2.0;
			if(random.nextDouble()>prob) sim.spins.set(i%sim.L, i/sim.L, -1);
			else sim.spins.set(i%sim.L, i/sim.L, 1);
		}

	}
	
	private double [] getSymmetricSlice(String file){
		double [] slice = new double [sim.L];
		double [] temp= FileUtil.readConfigFromFile(file, sim.L);
		double minPhi0Value = 1.0;
		int minPhi0Location = -1;
		for (int i = 0; i < sim.L; i++){
			if (temp[i] < minPhi0Value){
				minPhi0Location = i;
				minPhi0Value = temp[i];
			}
		}	
		for (int i = 0; i < sim.L; i++){
			slice[i] = temp[(minPhi0Location+i)%sim.L];
		}
		return slice;
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

	
	private void writeStSCtoFile(){
		String message1 = "#Glauber Monte Carlo run: S vs t for several values of k. Stripe to clump H quench.";
		String fileName = params.sget("Data Dir") + File.separator + "f00";
		StringBuffer fileBuffer = new StringBuffer(); fileBuffer.append(fileName);
		for (int i=0; i < accNo; i ++){
			StringBuffer mb = new StringBuffer();
			mb.append("# kx value = ");	mb.append(kRvalue(sfLabel[i]%sim.L)); mb.append(" ky value = ");
			mb.append(kRvalue(ky));			
			fileBuffer.deleteCharAt(fileBuffer.length()-1);fileBuffer.deleteCharAt(fileBuffer.length()-1);
			if(i < 10) fileBuffer.append("0");
			fileBuffer.append(i); fileName = fileBuffer.toString();
			String message2 = mb.toString();
			FileUtil.initFile(fileName, params, message1, message2);		
			FileUtil.printAccumToFile(fileName, sf_tAveAcc[i]);
		}
		
	}
	
	
	public void initialize(){
		ky = params.iget("ky");
		sim = new IsingLR(params);
		for (int i = 0; i < accNo; i++){
			sf_tAcc[i] = new Accumulator();
			sf_tAveAcc[i] = new Accumulator(); sf_tAveAcc[i].enableErrorBars(true);
			sfLabel[i] = ky*sim.L+i*params.iget("kx int");
			System.out.println("sfLabel " + i + " = " + sfLabel[i] + " kRValue hor = " + kRvalue(sfLabel[i]%sim.L));
		}
		sim.randomizeField(params.fget("Initial magnetization"));		
		dx = 1;
		String inputFile = params.sget("Input 1D File");
		inputSlice = getSymmetricSlice(inputFile);
		String unstableFile = params.sget("Input 1D Unstable File");
		unstableSlice = getSymmetricSlice(unstableFile);
	}
	
	public void writeConfigToFile(){
		String configFileName = "../../../research/javaData/stripeToClumpInvestigation/monteCarloData/configs/config";
		FileUtil.deleteFile(configFileName);
		FileUtil.writeConfigToFile(configFileName, (sim.L/dx)*(sim.L/dx), sim.getField(dx));
	}
	

		
}

package rachele.ising.dim2.apps;

import static java.lang.Math.PI;
import static java.lang.Math.abs;
import static scikit.util.Utilities.format;

import java.awt.Color;
import java.io.File;
import java.util.Random;

import rachele.ising.dim2.IsingLR;
import rachele.util.FileUtil;
import rachele.util.FourierTransformer;
import scikit.dataset.Accumulator;
import scikit.graphics.dim2.Grid;
import scikit.graphics.dim2.Plot;
import scikit.jobs.Control;
import scikit.jobs.Job;
import scikit.jobs.Simulation;
import scikit.jobs.params.ChoiceValue;
import scikit.jobs.params.DirectoryValue;
import scikit.jobs.params.FileValue;

/**
* 
* Monte Carlo Simulation to produce structure factor vs time averages.
* This should take an average of horizontal and vertical structure factor info.
* 
*/
public class MCStripesClumpsSt1DsolnApp extends Simulation{

	Grid grid = new Grid("Long Range Ising Model");
	Plot sftPlot = new Plot("sf_t plot");
	int dx;
	IsingLR sim;
	public FourierTransformer fft;
	double [] sFactor;
	Random random = new Random();
	Accumulator sf_kAcc;
	int accNo = 6;
	Accumulator [] sf_tAveAcc = new Accumulator [accNo];
	Accumulator [] sf_tAcc = new Accumulator [accNo];
	int [] sfLabel = new int [accNo];
	int ky;
	
	public static void main(String[] args) {
		new Control(new MCStripesClumpsSt1DsolnApp(), "Stripes to Clumps");
	}

	public void load(Control c) {

		c.frame(grid);
		c.frame(sftPlot);
		params.add("Data Dir",new DirectoryValue("/home/erdomi/data/lraim/stripeToClumpInvestigation/mcResults/conservedOP/testruns"));
		params.add("Input 1D File",new FileValue("/home/erdomi/data/lraim/configs1dAutoName/L128R45T0.08h0.65"));
		params.addm("Dynamics", new ChoiceValue("Kawasaki Glauber", "Ising Glauber", "Kawasaki Metropolis",  "Ising Metropolis"));
		params.add("Random seed", 0);
		params.add("L", 1<<7);
		params.add("R", 46);//1<<6);
		params.add("Initial magnetization", 0.0);
		params.addm("T", 0.04);
		params.addm("J", -1.0);
		params.addm("h", 0.0);
		params.addm("dt", 1.0);//1/(double)(1<<1));
		params.addm("maxTime", 30.0);
		params.addm("ky", 2);
		params.addm("kx int", 2);
		params.add("time");
		params.add("magnetization");
		params.add("Lp");
		params.add("Reps");
		flags.add("Write Config");
	}	

	public void animate() {
		sftPlot.setLogScale(false, true);
		params.set("time", format(sim.time()));
		params.set("magnetization", format(sim.magnetization()));
		sim.setParameters(params);
		params.set("Lp", sim.L/dx);
		grid.registerData(sim.L/dx, sim.L/dx, sim.getField(dx));
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
		int sfLabel = findBestkR();
		System.out.println(sfLabel);

		while (true) {
			for (int i = 0; i < accNo; i ++)
				sf_tAcc[i].clear();
			initializeStripes();
			sim.restartClock();
			sFactor = fft.calculate2DSF(sim.getField(dx), false, false);
			sim.restartClock();
			int recordInt = 0;
			int recordStep = 0;
			step = 0.25;
			while (sim.time() < maxTime){
				sim.step();
				Job.animate();
				if (sim.time() > recordStep){
					sFactor = fft.find2DSF(sim.getField(dx), sim.L);
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

	private void collect(double [] sFactor, int sfLabelHor){
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

		String fileName = params.sget("Input 1D File");
		//need to make phi0 symmetric
		double [] tempPhi0 = FileUtil.readConfigFromFile(fileName, sim.L);
		double [] phi0 = new double [sim.L];
		double minPhi0Value = 1.0;
		int minPhi0Location = -1;
		for (int i = 0; i < sim.L; i++){
			if (tempPhi0[i] < minPhi0Value){
				minPhi0Location = i;
				minPhi0Value = tempPhi0[i];
				//System.out.println(tempPhi0[i] + " " + i);
			}
		}	
		//System.out.println(tempPhi0[minPhi0Location] + " " + minPhi0Location);
		for (int i = 0; i < sim.L; i++){
			phi0[i] = tempPhi0[(minPhi0Location+i)%sim.L];
			//System.out.println("phi0 " + i + " = " + phi0[i]);
		}		
		//now load the Ising lattice with the proper probability
		for (int i = 0; i < sim.L*sim.L; i++){
			double prob = (phi0[i%sim.L]+1.0)/2.0;
			if(random.nextDouble()>prob) sim.spins.set(i%sim.L, i/sim.L, -1);
			else sim.spins.set(i%sim.L, i/sim.L, 1);
		}

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
		String fileName = params.sget("Data Dir") + File.separator + "f0";
		StringBuffer fileBuffer = new StringBuffer(); fileBuffer.append(fileName);
		for (int i=0; i < accNo; i ++){
			StringBuffer mb = new StringBuffer();
			mb.append("# kx value = ");	mb.append(kRvalue(sfLabel[i]%sim.L)); mb.append(" ky value = ");
			mb.append(kRvalue(ky));			
			fileBuffer.deleteCharAt(fileBuffer.length()-1);	fileBuffer.append(i); fileName = fileBuffer.toString();
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
	}
	
	public void writeConfigToFile(){
		String configFileName = "../../../research/javaData/stripeToClumpInvestigation/monteCarloData/configs/config";
		FileUtil.deleteFile(configFileName);
		FileUtil.writeConfigToFile(configFileName, (sim.L/dx)*(sim.L/dx), sim.getField(dx));
	}
	

		
}

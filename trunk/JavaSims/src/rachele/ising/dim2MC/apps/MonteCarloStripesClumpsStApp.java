package rachele.ising.dim2MC.apps;

import static java.lang.Math.PI;
import static java.lang.Math.abs;
import static scikit.util.Utilities.format;
import java.awt.Color;
import java.io.DataInputStream;
import java.io.EOFException;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.util.Random;

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
import scikit.jobs.params.DirectoryValue;
import scikit.jobs.params.FileValue;

/**
* 
* Monte Carlo Simulation to produce structure factor vs time averages.
* This has code to characterize if the stripes are hor or vertical and
* collects only the dominant direction data.
* 
* For 1D input, use MCStripesClumpsSt1DsolnApp
* 
*/

public class MonteCarloStripesClumpsStApp extends Simulation{

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
	boolean init1D;	
	
	public static void main(String[] args) {
		new Control(new MonteCarloStripesClumpsStApp(), "Monte Carlo");
	}

	public void load(Control c) {

		c.frame(grid);
		c.frame(sftPlot);
		params.add("Data Dir",new DirectoryValue("/Users/erdomi/data/lraim/stripeToClumpInvestigation/mcResults/noConservedOP/r100"));
		params.add("Input 1D File",new FileValue("/Users/erdomi/data/lraim/configs1D/L128R50T0-04h0"));
		params.addm("Dynamics", new ChoiceValue("Ising Glauber","Kawasaki Glauber", "Kawasaki Metropolis",  "Ising Metropolis"));
		params.addm("init", new ChoiceValue( "Use 1D Soln", "Init Stripes from random"));
		params.add("Random seed", 0);
		params.add("L", 1<<7);
		params.add("R", 46);//1<<6);
		params.addm("Jump Range", 46);
		params.add("Initial magnetization", 0.0);
		params.addm("T", 0.04);
		params.addm("J", -1.0);
		params.addm("h", 0.8);
		params.addm("dt", 1/(double)(1<<3));
		params.addm("maxTime", 10.0);
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
		double initTime = 15.0;
		double maxTime = params.fget("maxTime");//max time after quench time
		double quenchH = params.fget("h");// ext field after quench (initialized at h=0)
		int repNo = 0;
		int sfLabel = findBestkR();
		System.out.println(sfLabel);

		while (true) {
			for (int i = 0; i < accNo; i ++)
				sf_tAcc[i].clear();
			sim.restartClock();
			initializeStripes(initTime, quenchH);
			sim.restartClock();
			sFactor = fft.calculate2DSF(sim.getField(dx), false, false);
			boolean vertStripes;
			if(sFactor[sfLabel]>sFactor[sfLabel*sim.L/dx])vertStripes = true;
			else vertStripes = false;
			sim.restartClock();
			int recordInt = 0;
			int recordStep = 0;
			step = 0.25;
			while (sim.time() < maxTime){
				sim.step();
				Job.animate();
				if (sim.time() > recordStep){
					sFactor = fft.find2DSF(sim.getField(dx), sim.L);
					collect(sFactor, vertStripes,sfLabel);
					recordStep += step;
					recordInt +=1;
				}
			}	
			repNo += 1;
			params.set("Reps", repNo);
			writeStSCtoFile(sfLabel, initTime, kRvalues(vertStripes, sfLabel));
		}
	}

	private void collect(double [] sFactor, boolean vertStripes, int sfLabelHor){
		int sfLabelVert = sfLabelHor*sim.L/dx;
		if(vertStripes){
			for (int i = 0; i < accNo; i ++)
				sfLabel[i] = sfLabelVert + i;
		}else{
			for (int i = 0; i < accNo; i ++)
				sfLabel[i] = sfLabelHor + i*sim.L/dx;
		}
		for (int i = 0; i < accNo; i ++){
			//int kx = sfLabel[i]%sim.L/dx; int ky = sfLabel[i]/sim.L/dx;
			//double kxR = 2*PI*kx*sim.R/sim.L;double kyR = 2*PI*ky*sim.R/sim.L;
			//System.out.println("kxR = " + kxR + " kyR =  " + kyR);
			sf_tAveAcc[i].accum(sim.time(),sFactor[sfLabel[i]]);			
			sf_tAcc[i].accum(sim.time(),sFactor[sfLabel[i]]);
		}
	}
	
	private double [] kRvalues(boolean vertStripes, int sfLabelHor){
		double [] kRvalue = new double [accNo*2];
		int [] kRCoord = kRLabels(vertStripes, sfLabelHor);
		for (int i = 0; i < accNo; i++){
			kRvalue[2*i] = 2*PI*kRCoord[2*i]*sim.R/(sim.L);
			kRvalue[2*i+1] = 2*PI*kRCoord[2*i+1]*sim.R/(sim.L); 
		}	
		return kRvalue;
	}
		
		private int [] kRLabels(boolean vertStripes, int sfLabelHor){
			int [] kRCoord = new int [2*accNo];
			int sfLabelVert = sfLabelHor*sim.L/dx;
			if(vertStripes){
				for (int i = 0; i < accNo; i ++){
					int label = sfLabelVert + i;
					kRCoord[i*2] = label %sim.L/dx;
					kRCoord[i*2+1] = label /(sim.L/dx);
				}
			}else{
				for (int i = 0; i < accNo; i ++){
					int label = sfLabelHor + i*sim.L/dx;
					kRCoord[i*2] = label %sim.L/dx;
					kRCoord[i*2+1] = label /(sim.L/dx);					
				}
			}
			return kRCoord;
		}
	
	private void initializeStripes(double initTime, double finalH){
		params.set("h", 0.0);
		sim.randomizeField(params.fget("Initial magnetization"));	
		while(sim.time() < initTime){ 
			sim.step();
			Job.animate();		
		}
		params.set("h", finalH);
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

	
	private void writeStSCtoFile(int sfInt, double initializeTime, double [] kRvalues){
		String message1 = "#Glauber Monte Carlo run: S vs t for several values of k. Stripe to clump H quench.";
		String fileName = params.sget("Data Dir") + File.separator + "f0";
		StringBuffer fileBuffer = new StringBuffer(); fileBuffer.append(fileName);
		for (int i=0; i < accNo; i ++){
			//System.out.println("start " + i);
			StringBuffer mb = new StringBuffer();
			//mb.append("# init time = "); mb.append(initializeTime);
			mb.append("# kx value = ");	mb.append(kRvalues[2*i]); mb.append(" ky value = ");
			//double krvalue = 2*sim.R*Math.PI*(sfInt+i)/sim.L;
			mb.append(kRvalues[2*i+1]);			
			fileBuffer.deleteCharAt(fileBuffer.length()-1);	fileBuffer.append(i); fileName = fileBuffer.toString();
			String message2 = mb.toString();
			FileUtil.initFile(fileName, params, message1, message2);		
			FileUtil.printAccumToFile(fileName, sf_tAveAcc[i]);
		}
		
	}
	
	
	public void initialize(){
		for (int i = 0; i < accNo; i++){
			sf_tAcc[i] = new Accumulator();
			sf_tAveAcc[i] = new Accumulator(); sf_tAveAcc[i].enableErrorBars(true);
		}
		sim = new IsingLR(params);
		sim.randomizeField(params.fget("Initial magnetization"));		
		dx = 1;
		if(params.sget("init") == "Read From File") readInitialConfiguration();
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
		
}

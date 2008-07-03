package rachele.ising.dim2.apps;

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


//Monte Carlo Structure Factor vs time with 1D initialization
public class MonteCarloStripesClumpsStApp extends Simulation{

	Grid grid = new Grid("Long Range Ising Model");
	Plot sftPlot = new Plot("sf_t plot");
	int dx;
	IsingLR sim;
	public FourierTransformer fft;
	double [] sFactor;
	Random random = new Random();
	Accumulator sf_kAcc;
	Accumulator [] sf_tAveAcc;
	Accumulator [] sf_tAcc;
	int [] sfLabel;
	int accNo;
	Accumulator sf_kTheoryAcc, sf_kTheory2Acc, sf_tTheoryAcc, sf_tTheory2Acc;

	public static void main(String[] args) {
		new Control(new MonteCarloStripesClumpsStApp(), "Monte Carlo");
	}

	public void load(Control c) {

		c.frame(grid);
		c.frame(sftPlot);
		params.addm("Dynamics", new ChoiceValue("Ising Glauber","Kawasaki Glauber", "Kawasaki Metropolis",  "Ising Metropolis"));
		params.addm("init", new ChoiceValue( "Random", "Read From File"));
		params.add("Random seed", 0);
		params.add("L", 1<<7);
		params.add("R", 10);//1<<6);
		params.add("Initial magnetization", 0.0);
		params.addm("T", 0.1);
		params.addm("J", -1.0);
		params.addm("h", 0.0);
		params.addm("dt", 1/(double)(1<<3));
		params.addm("maxTime", 5.0);
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
			sftPlot.registerLines(sb.toString(), sf_tAveAcc[i], Color.BLACK);//getColor("fill", i*2));
			sb = new StringBuffer();sb.append("s(t) "); sb.append(i);
			sftPlot.registerLines(sb.toString(), sf_tAcc[i], Color.getColor("fill", i*2+1));
		}
		sftPlot.registerLines("Blah", sf_tAcc[0], Color.BLACK);
		
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
		double initTime = 15.0;
		boolean init1D = true;
		double maxTime = 15.0;//max time after quench time
		double quenchH = 0.9;
		int repNo = 0;
		int sfLabel = findBestkR();

		while (true) {
			for (int i = 0; i < accNo; i ++)
				sf_tAcc[i].clear();
			initializeStripes(init1D, initTime);
			sim.randomizeField(0.0);
			sim.restartClock();
			sFactor = fft.calculate2DSF(sim.getField(dx), false, false);
			boolean vertStripes;
			if(sFactor[sfLabel]>sFactor[sfLabel*sim.L/dx])vertStripes = true;
			else vertStripes = false;
			sim.restartClock();
			params.set("h", quenchH);
			int recordInt = 0;
			int recordStep = 0;
			step = 0.25;
			while (sim.time() < maxTime){
				sim.step();
				Job.animate();
				if (sim.time() > recordStep){
					sFactor = fft.calculate2DSF(sim.getField(dx), false, false);
					collect(sFactor, vertStripes,sfLabel);
					recordStep += step;
					recordInt +=1;
				}
			}	
			repNo += 1;
			params.set("Reps", repNo);
			writeStSCtoFile(sfLabel, initTime);
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
			sf_tAveAcc[i].accum(sim.time(),sFactor[sfLabel[i]]);			
			sf_tAcc[i].accum(sim.time(),sFactor[sfLabel[i]]);
		}
	}
	
	private void initializeStripes(boolean init1D, double initTime){
		if(init1D){
			String fileName = "../../../research/javaData/configs1d/config";
			//need to make phi0 symmetric
			double [] tempPhi0 = FileUtil.readConfigFromFile(fileName, sim.L);
			double [] phi0 = new double [sim.L];
			double minPhi0Value = 1.0;
			int minPhi0Location = -1;
			for (int i = 0; i < sim.L; i++){
				if (tempPhi0[i] < minPhi0Value){
					minPhi0Location = i;
					minPhi0Value = tempPhi0[i];
					System.out.println(tempPhi0[i] + " " + i);
				}
			}	
			System.out.println(tempPhi0[minPhi0Location] + " " + minPhi0Location);
			for (int i = 0; i < sim.L; i++){
				phi0[i] = tempPhi0[(minPhi0Location+i)%sim.L];
				System.out.println("phi0 " + i + " = " + phi0[i]);
			}		
			//now load the Ising lattice with the proper probability
			for (int i = 0; i < sim.L*sim.L; i++){
				double prob = (phi0[i%sim.L]+1.0)/2.0;
				if(random.nextDouble()>prob) sim.spins.set(i%sim.L, i/sim.L, +1);
				else sim.spins.set(i%sim.L, i/sim.L, -1);
			}
		}else{
			params.set("h", 0.0);
			while(sim.time() < initTime){ 
				sim.step();
				Job.animate();
			}			
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

	
	private void writeStSCtoFile(int sfInt, double initializeTime){
		String message1 = "#Glauber Monte Carlo run: S vs t for several values of k. Stripe to clump H quench. Init H = 0.";
		String fileName = "../../../research/javaData/stripeToClumpInvestigation/monteCarloData/squareResults/svtSCinit1D1/s0";
		StringBuffer fileBuffer = new StringBuffer(); fileBuffer.append(fileName);
		for (int i=0; i < accNo; i ++){
			StringBuffer mb = new StringBuffer();
			mb.append("# init time = "); mb.append(initializeTime); mb.append(", sf label = ");	mb.append(sfInt+i); mb.append(" kR value = ");
			double krvalue = 2*sim.R*Math.PI*(sfInt+i)/sim.L;
			mb.append(krvalue);			
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

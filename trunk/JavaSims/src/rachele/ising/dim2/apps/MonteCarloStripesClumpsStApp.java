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
	Accumulator sf_tAveAcc0,sf_tAcc0;
	Accumulator sf_tAveAcc1,sf_tAcc1;
	Accumulator sf_tAveAcc2,sf_tAcc2;
	Accumulator sf_tAveAcc3,sf_tAcc3;
	Accumulator sf_tAveAcc4,sf_tAcc4;
	Accumulator sf_tAveAcc5,sf_tAcc5;
	Accumulator sf_kTheoryAcc, sf_kTheory2Acc, sf_tTheoryAcc, sf_tTheory2Acc;
	double [] sfTimeArray;

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

		sftPlot.registerLines("st(k) Ave0", sf_tAveAcc0, Color.BLACK);
		sftPlot.registerLines("st(k)0", sf_tAcc0, Color.BLACK);

		sftPlot.registerLines("st(k) Ave1", sf_tAveAcc1, Color.LIGHT_GRAY);
		sftPlot.registerLines("st(k)1", sf_tAcc1, Color.LIGHT_GRAY);

		sftPlot.registerLines("st(k) Ave2", sf_tAveAcc2, Color.BLUE);
		sftPlot.registerLines("st(k)2", sf_tAcc2, Color.BLUE);

		sftPlot.registerLines("st(k) Ave3", sf_tAveAcc3, Color.CYAN);
		sftPlot.registerLines("st(k)3", sf_tAcc3, Color.CYAN);

		sftPlot.registerLines("st(k) Ave4", sf_tAveAcc4, Color.GREEN);
		sftPlot.registerLines("st(k)4", sf_tAcc4, Color.GREEN);

		sftPlot.registerLines("st(k) Ave5", sf_tAveAcc5, Color.RED);
		sftPlot.registerLines("st(k)5", sf_tAcc5, Color.RED);

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
		sf_tAcc0 = new Accumulator();
		sf_tAcc1 = new Accumulator();
		sf_tAcc2 = new Accumulator();
		sf_tAcc3 = new Accumulator();
		sf_tAcc4 = new Accumulator();
		sf_tAcc5 = new Accumulator();
		sf_tAveAcc0 = new Accumulator();
		sf_tAveAcc1 = new Accumulator();
		sf_tAveAcc2 = new Accumulator();
		sf_tAveAcc3 = new Accumulator();
		sf_tAveAcc4 = new Accumulator();
		sf_tAveAcc5 = new Accumulator();

		
		int sfLabel = findBestkR();

		sfTimeArray = new double[(int)((maxTime/step)+1)];

		while (true) {
			sf_tAcc0.clear();
			sf_tAcc1.clear();
			sf_tAcc2.clear();
			sf_tAcc3.clear();
			sf_tAcc4.clear();
			sf_tAcc5.clear();
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
		int sfLabel0, sfLabel1, sfLabel2, sfLabel3, sfLabel4, sfLabel5;
		int sfLabelVert = sfLabelHor*sim.L/dx;
		if(vertStripes){
			sfLabel0 = sfLabelVert;
			sfLabel1 = sfLabelVert + 1;
			sfLabel2 = sfLabelVert + 2;
			sfLabel3 = sfLabelVert + 3;
			sfLabel4 = sfLabelVert + 4;
			sfLabel5 = sfLabelVert + 5;
		}else{
			sfLabel0 = sfLabelHor;
			sfLabel1 = sfLabelHor + 1*sim.L/dx;
			sfLabel2 = sfLabelHor + 2*sim.L/dx;
			sfLabel3 = sfLabelHor + 3*sim.L/dx;
			sfLabel4 = sfLabelHor + 4*sim.L/dx;
			sfLabel5 = sfLabelHor + 5*sim.L/dx;
		}
		sf_tAveAcc0.accum(sim.time(),sFactor[sfLabel0]);
		sf_tAveAcc1.accum(sim.time(),sFactor[sfLabel1]);
		sf_tAveAcc2.accum(sim.time(),sFactor[sfLabel2]);
		sf_tAveAcc3.accum(sim.time(),sFactor[sfLabel3]);
		sf_tAveAcc4.accum(sim.time(),sFactor[sfLabel4]);
		sf_tAveAcc5.accum(sim.time(),sFactor[sfLabel5]);
		sf_tAcc0.accum(sim.time(),sFactor[sfLabel0]);
		sf_tAcc1.accum(sim.time(),sFactor[sfLabel1]);
		sf_tAcc2.accum(sim.time(),sFactor[sfLabel2]);
		sf_tAcc3.accum(sim.time(),sFactor[sfLabel3]);
		sf_tAcc4.accum(sim.time(),sFactor[sfLabel4]);
		sf_tAcc5.accum(sim.time(),sFactor[sfLabel5]);		
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

		StringBuffer messageBuffer = new StringBuffer();
		messageBuffer.append("# sf label = ");	messageBuffer.append(sfInt); messageBuffer.append(" kR values = ");
		double krvalue = 2*sim.R*Math.PI*sfInt/sim.L;
		messageBuffer.append(krvalue);	messageBuffer.append(", ");
		krvalue = 2*sim.R*Math.PI*(sfInt+1)/sim.L;
		messageBuffer.append(krvalue);	messageBuffer.append(", ");
		krvalue = 2*sim.R*Math.PI*(sfInt+2)/sim.L;
		messageBuffer.append(krvalue);	messageBuffer.append(", ");
		krvalue = 2*sim.R*Math.PI*(sfInt+3)/sim.L;
		messageBuffer.append(krvalue);	messageBuffer.append(", ");
		krvalue = 2*sim.R*Math.PI*(sfInt+4)/sim.L;
		messageBuffer.append(krvalue);	messageBuffer.append(", ");
		krvalue = 2*sim.R*Math.PI*(sfInt+5)/sim.L;
		messageBuffer.append(krvalue);messageBuffer.append(", init time =");
		messageBuffer.append(initializeTime);
		String message1 = "#Glauber Monte Carlo run: S vs t for several values of k. Stripe to clump H quench. Init H = 0.";
		String message2 = messageBuffer.toString();
		String fileName = "../../../research/javaData/stripeToClumpInvestigation/monteCarloData/squareResults/svtRunQuenchH2/s";
		FileUtil.initFile(fileName, params, message1, message2);
		FileUtil.printAccumToFile(fileName, sf_tAveAcc0);
		
		StringBuffer fileBuffer = new StringBuffer(); fileBuffer.append(fileName);fileBuffer.append(1);
		String file1 = fileBuffer.toString();
		FileUtil.initFile(file1, params, message1, message2);		
		FileUtil.printAccumToFile(file1, sf_tAveAcc1);
		
		fileBuffer.deleteCharAt(fileBuffer.length()-1);	fileBuffer.append(2);String file2 = fileBuffer.toString();
		FileUtil.initFile(file2, params, message1, message2);		
		FileUtil.printAccumToFile(file2, sf_tAveAcc2);

		fileBuffer.deleteCharAt(fileBuffer.length()-1);	fileBuffer.append(3);String file3 = fileBuffer.toString();
		FileUtil.initFile(file3, params, message1, message2);		
		FileUtil.printAccumToFile(file3, sf_tAveAcc3);
		
		fileBuffer.deleteCharAt(fileBuffer.length()-1);	fileBuffer.append(4);String file4 = fileBuffer.toString();
		FileUtil.initFile(file4, params, message1, message2);		
		FileUtil.printAccumToFile(file4, sf_tAveAcc4);
		
		fileBuffer.deleteCharAt(fileBuffer.length()-1);	fileBuffer.append(5);String file5 = fileBuffer.toString();
		FileUtil.initFile(file5, params, message1, message2);		
		FileUtil.printAccumToFile(file5, sf_tAveAcc5);
	}
	
	
	public void initialize(){
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

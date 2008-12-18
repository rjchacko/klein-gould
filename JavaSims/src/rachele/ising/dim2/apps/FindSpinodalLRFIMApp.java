package rachele.ising.dim2.apps;

import static scikit.util.Utilities.format;

import java.awt.Color;
import java.io.File;

import rachele.ising.dim2.IsingLR;
import rachele.util.FileUtil;
import scikit.dataset.Accumulator;
import scikit.graphics.dim2.Grid;
import scikit.graphics.dim2.Plot;
import scikit.jobs.Control;
import scikit.jobs.Job;
import scikit.jobs.Simulation;
import scikit.jobs.params.ChoiceValue;
import scikit.jobs.params.DirectoryValue;
/**
 * Finds effective spinodal for the long-range Ising model with Glauber dynamics using short time universal scaling.
 */
public class FindSpinodalLRFIMApp extends Simulation{

	Grid grid = new Grid("Long Range Ising Model");
	Plot magPlot = new Plot("Magnetization");
	Plot magSqPlot = new Plot("Magnetization2");
	Plot chiPlot = new Plot("Susceptibility");
	Plot chiTPlot = new Plot("Susceptibility v T");
	IsingLR sim;
//	IsingArbitraryConnect sim;
	boolean measure=false;
	int dx,n, repNo;
	double mag, maxTime, delta_mag;
	double [] aveMag, aveMagSq; 
	Accumulator magAveAcc = new Accumulator();
	Accumulator magSqAveAcc = new Accumulator();
	Accumulator chiAveAcc = new Accumulator();
	Accumulator magAcc = new Accumulator();
	Accumulator magSqAcc = new Accumulator();

	public static void main(String[] args) {
		new Control(new FindSpinodalLRFIMApp(), "Monte Carlo");
	}

	public void load(Control c) {
		c.frameTogether("MC data",grid,magPlot,magSqPlot);
		params.add("Data Dir",new DirectoryValue("/home/erdomi/data/spinodal_find/testRuns"));
		params.addm("Dynamics", new ChoiceValue("Ising Glauber","Kawasaki Glauber", "Kawasaki Metropolis",  "Ising Metropolis"));
		params.add("Random seed", 0);
		params.add("L", 1<<7);
		params.add("R", 6);//1<<6);
		params.add("Initial magnetization", 1.0);
		params.addm("T", 0.4444444444);
		params.addm("J", 1.0);
		params.addm("h", -0.31);
		params.addm("dt", 1.0);//1/(double)(1<<4));
		params.addm("take data",1);
		params.addm("max time",100);
		params.add("rep no");
		params.add("time");
		params.add("magnetization");
		flags.add("Clear Accs");
	}


	public void animate() {
		grid.setScale(-1.0, 1.0);
		grid.registerData(sim.L, sim.L, sim.getField(1));
		magPlot.setLogScale(true, true);
		magSqPlot.setLogScale(true, true);
		maxTime = params.fget("max time");
		params.set("time", format(sim.time()));
		params.set("magnetization", format(mag));
		params.set("rep no", repNo);
		
		magPlot.registerLines("Magnetization", magAveAcc, Color.black);
		magPlot.registerLines("Magnetization1", magAcc, Color.green);
		magSqPlot.registerLines("secMom", magSqAveAcc, Color.black);
		magSqPlot.registerLines("secMom1", magSqAcc, Color.green);


		if(flags.contains("Clear Accs")){
			magAveAcc.clear();
			magSqAveAcc.clear();
			chiAveAcc.clear();
		}
		if(flags.contains("Measure")){
			if (measure) measure=false;
			else if (measure==false) measure=true;
		}
		flags.clear();
		sim.setParameters(params);
	}


	public void clear() {


	}


	public void run() {
		magAveAcc.enableErrorBars(true);
		magSqAveAcc.enableErrorBars(true);
		magAveAcc.clear();
		magSqAveAcc.clear();
		chiAveAcc.clear();
		sim = new IsingLR(params);
		maxTime = params.fget("max time");
		repNo = 0;
		double mag_sp = Math.sqrt(1.0-sim.T/sim.J);
		
		int maxTimeChunk = (int)(maxTime/sim.dTime());
		aveMag = new double [maxTimeChunk];
		aveMagSq = new double [maxTimeChunk];
//		sim = new IsingArbitraryConnect(params);
		
		while(true){
			sim.restartClock();
			magAcc.clear();
			sim.randomizeField(params.fget("Initial magnetization"));
			repNo += 1;
			mag = sim.magnetization();
			for(int i = 0; i < maxTimeChunk; i ++){
				if(mag > 0){
					sim.step();
					mag = sim.magnetization();
					delta_mag = mag - mag_sp;
					magAveAcc.accum(sim.time(),delta_mag);
					magAcc.accum(sim.time(),delta_mag);
					aveMag[i] = (aveMag[i]*(repNo-1)+delta_mag)/(repNo);
					aveMagSq[i] = (aveMagSq[i]*(repNo-1)+delta_mag*delta_mag)/(repNo);
					magSqAveAcc.accum(sim.time(),aveMagSq[i]-aveMag[i]*aveMag[i]);
					Job.animate();
				}
			}
			Job.animate();
			writeToFile();
		}


	}

	void writeToFile(){	
		String message1 = "#Glauber Monte Carlo run: Short time scaling dynamics to find spinodal field h_s and critical exponents.";
		String fileNameM = params.sget("Data Dir") + File.separator + "h" + -sim.h +"M";
		StringBuffer fileBuffer = new StringBuffer(); fileBuffer.append(fileNameM);
		StringBuffer mb = new StringBuffer(); 
		mb.append("#Average magnetization vs time");
		String message2 = mb.toString();
		FileUtil.initFile(fileNameM, params, message1, message2);		
		FileUtil.printAccumToFile(fileNameM, magAveAcc);
		
		String fileNameM2 = params.sget("Data Dir") + File.separator + "h" + -sim.h +"M2";		
		StringBuffer mb2 = new StringBuffer();
		mb2.append("#2nd moment of mag vs time");
		String message22 = mb2.toString();
		FileUtil.initFile(fileNameM2, params, message1, message22);		
		FileUtil.printAccumToFile(fileNameM2, magSqAveAcc);
	}
	
}



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

public class MCspinodalFind extends Simulation{

	Grid grid = new Grid("Long Range Ising Model");
	Plot magPlot = new Plot("Magnetization");
	Plot magSqPlot = new Plot("Magnetization2");
	Plot chiPlot = new Plot("Susceptibility");
	Plot chiTPlot = new Plot("Susceptibility v T");
	IsingLR sim;
//	IsingArbitraryConnect sim;
	boolean measure=false;
	int dx,n, repNo;
	double mag, maxTime;
	double [] aveMag, aveMagSq; 
	Accumulator magAveAcc = new Accumulator();
	Accumulator magSqAveAcc = new Accumulator();
	Accumulator chiAveAcc = new Accumulator();
	Accumulator magAcc = new Accumulator();
	Accumulator magSqAcc = new Accumulator();
	Accumulator chiAcc = new Accumulator();
	Accumulator chiTAcc = new Accumulator();

	public static void main(String[] args) {
		new Control(new MCspinodalFind(), "Monte Carlo");
	}

	public void load(Control c) {
		c.frameTogether("MC data",grid,magPlot,magSqPlot);
		params.add("Data Dir",new DirectoryValue("/home/erdomi/data/spinodal_find/testRuns"));
		params.addm("Dynamics", new ChoiceValue("Ising Glauber","Kawasaki Glauber", "Kawasaki Metropolis",  "Ising Metropolis"));
		params.add("Random seed", 0);
		params.add("L", 1<<9);
		params.add("R", 184);//1<<6);
		params.add("Initial magnetization", 1.0);
		params.addm("T", 0.4444444444);
		params.addm("J", 1.0);
		params.addm("h", -0.317612);
		params.addm("dt", 0.1);//1/(double)(1<<4));
		params.addm("take data",1);
		params.addm("max time",30);
		params.add("rep no");
		params.add("time");
		params.add("magnetization");
		flags.add("Clear Accs");
	}


	public void animate() {
		grid.setScale(-1.0, 1.0);
		grid.registerData(sim.L/dx, sim.L/dx, sim.getField(dx));
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
//		chiPlot.registerLines("Susceptibility", chiAveAcc, Color.blue);
//		chiPlot.registerLines("Susceptibility1", chiAcc, Color.green);

		if(flags.contains("Clear Accs")){
			magAveAcc.clear();
			magSqAveAcc.clear();
			chiAveAcc.clear();
		}
		if(flags.contains("Measure")){
			if (measure) measure=false;
			else if (measure==false) measure=true;
		}
//		int measureStep = params.iget("take data");
//		params.set("Measuring", measure);
//		if(measure){
//			if(sim.time()%measureStep==0){
//				chiTAcc.accum(sim.T,chi);
//				chiTPlot.registerLines("chi T", chiTAcc, Color.red);
//			}
//		}
		flags.clear();
		sim.setParameters(params);
	}


	public void clear() {


	}


	public void run() {
		chiAcc.enableErrorBars(true);
		magAveAcc.enableErrorBars(true);
		magSqAveAcc.enableErrorBars(true);
		magAveAcc.clear();
		magSqAveAcc.clear();
		chiAveAcc.clear();
		sim = new IsingLR(params);
		maxTime = params.fget("max time");
		repNo = 0;
		
		int maxTimeChunk = (int)(maxTime/sim.dTime());
		aveMag = new double [maxTimeChunk];
		aveMagSq = new double [maxTimeChunk];
//		sim = new IsingArbitraryConnect(params);
		dx=1;
		while(true){
			sim.restartClock();
			magSqAcc.clear();
			chiAcc.clear();
			magAcc.clear();
			sim.randomizeField(params.fget("Initial magnetization"));
			repNo += 1;
			for(int i = 0; i < maxTimeChunk; i ++){
				sim.step();
				mag = sim.magnetization();
				magAveAcc.accum(sim.time(),mag);
				magAcc.accum(sim.time(),mag);
				aveMag[i] = (aveMag[i]*n+mag)/(n+1);
				aveMagSq[i] = (aveMagSq[i]*n+mag*mag)/(n+1);
				n+=1;
//				double chi = (aveMagSq[i]-aveMag[i]*aveMag[i])/sim.T;
//				chiAcc.accum(sim.time(), chi);
//				chiAveAcc.accum(sim.time(), chi);
				double secMom = (aveMagSq[i]-aveMag[i]*aveMag[i]);
				magSqAveAcc.accum(sim.time(),secMom);
				magSqAcc.accum(sim.time(),secMom);				
				Job.animate();
			}
			writeToFile();
		}

	}

	void writeToFile(){	
		String message1 = "#Glauber Monte Carlo run: Short time scaling dynamics to find spinodal field h_s.";
		String fileName = params.sget("Data Dir") + File.separator + "h" + -sim.h +"M";
		StringBuffer fileBuffer = new StringBuffer(); fileBuffer.append(fileName);
		StringBuffer mb = new StringBuffer(); String message2 = mb.toString();
		mb.append("Average magnetization vs time");
		FileUtil.initFile(fileName, params, message1, message2);		
		FileUtil.printAccumToFile(fileName, magAveAcc);
	}
	
}



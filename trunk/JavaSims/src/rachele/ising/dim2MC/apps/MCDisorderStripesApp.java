package rachele.ising.dim2MC.apps;

import static java.lang.Math.PI;
import static scikit.util.Utilities.format;
import java.awt.Color;
import java.io.File;

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

/**
* 
* Monte Carlo Simulation to produce structure factor vs time averages.
* Separates data into dominant and non-dominant directions for several k values
* to investigate disorder -> clump -> stripe transition.
*/


public class MCDisorderStripesApp extends Simulation{
	Grid grid = new Grid("Long Range Ising Model");
	Grid sfGrid = new Grid("SF");
	Plot sfkPlot = new Plot("sf_k plot");
	Plot sftPlot = new Plot("sf_t plot");
	Plot absMag = new Plot("abs mag");
	int dx; 
	IsingLR sim;
	public FourierTransformer fft;
	double [] sFactor;
	Accumulator sf_kAcc; 

	Accumulator  sf_tDomAveAcc =  new Accumulator(); //Average sf for dominant direction (H or V)
	Accumulator  sf_tNDomAveAcc = new Accumulator();//Average sf for non-dominant direction (H or V)
	Accumulator mag = new Accumulator();
	double [][] sf_tAcc;// 1st label: 0 = Horizontal direction, 1 = vertical direction; 2nd label -> time Label
    boolean clearFile;
   	double [] sfTimeArray;
	
	public static void main(String[] args) {
		new Control(new MCDisorderStripesApp(), "MC Disorder -> Stripes");
	}
	
	public void load(Control c) {
		
		c.frame(grid);
		c.frameTogether("plots", sftPlot, absMag);
		params.add("Data Dir",new DirectoryValue("/Users/erdomi/data/lraim/stripeToClumpInvestigation/mcResults/noConservedOP/testRuns"));
		params.addm("Dynamics", new ChoiceValue("Ising Glauber","Kawasaki Glauber", "Kawasaki Metropolis",  "Ising Metropolis"));
		params.add("Random seed", 0);
		params.add("L", 1<<8);
		params.add("R", 92);//1<<6);
		params.add("Initial magnetization", 0.0);
		params.addm("T", 0.096548444);
		params.addm("J", -1.0);
		params.addm("h", 0.0);
		params.addm("ky", 0);
		params.addm("dkx", 0);
		params.addm("dt", 1/(double)(1<<5));
		params.addm("maxTime", 20);
		params.add("time");
		params.add("magnetization");
		params.add("Lp");
		params.add("Reps");
	}	
	
	public void animate() {
		sftPlot.setLogScale(false, true);
		params.set("time", format(sim.time()));
		params.set("magnetization", format(sim.magnetization()));
		sim.setParameters(params);
		params.set("Lp", sim.L/dx);
		grid.registerData(sim.L/dx, sim.L/dx, sim.getField(dx));
		sfGrid.registerData(sim.L/dx, sim.L/dx, sFactor);
			StringBuffer sb1 = new StringBuffer();sb1.append("s(t) Dominant Direction"); 
			StringBuffer sb2 = new StringBuffer();sb2.append("s(t) NonDominant Direction");
			sftPlot.registerLines(sb1.toString(), sf_tDomAveAcc, Color.red);
			sftPlot.registerLines(sb2.toString(), sf_tNDomAveAcc, Color.blue);
		absMag.registerLines("ab mag", mag, Color.black);
	}
	
	public void clear() {
	}
	
	public void run() {
		initialize();
		mag.clear();
		fft = new FourierTransformer((int)(sim.L/dx));
		sFactor = new double [sim.L/dx*sim.L/dx];
		double maxTime = params.fget("maxTime");
		double recordStep = 0.25;
		double mod = maxTime%recordStep;
		int maxTimeLabel = (int)(maxTime/recordStep);
		System.out.println("Record Step = " + recordStep + " mod (must be zero) = " + mod);
		sf_tAcc = new double [2][maxTimeLabel];
		int repNo = 0;

		while (true) {
			//must clear all arrays because we will add 4 contributions to average.
			for (int i = 0; i < 2; i++){
					for (int k = 0; k < maxTimeLabel; k ++){
						sf_tAcc[i][k] = 0.0;						
					}
			}

			sim.randomizeField(params.fget("Initial magnetization"));
			sim.restartClock();
			int timeLabel =0;
			while (timeLabel < maxTimeLabel){
				Job.animate();
				if (sim.time()% recordStep == 0){
					sFactor = fft.calculate2DSF(sim.getField(dx), false, false);
					mag.accum(sim.time(), Math.abs(sim.magnetization()));
					sf_tAcc[0][timeLabel] = sFactor[2];//horizontal direction
					sf_tAcc[1][timeLabel] = sFactor[sim.L*2];//vertical direction
					timeLabel += 1;
				}
				sim.step();
			}	
			// find hor or vert stripes
			int domInt,nonDomInt;
			if(sFactor[sim.L*2]>sFactor[2]){
				domInt = 1; nonDomInt = 0;
				System.out.println("horizontal stripes");
			}else{
				domInt = 0; nonDomInt = 1;
				System.out.println("vertical stripes");
			}
			// record in dom and non-dominant
//			for (int i = 0; i < accNo; i++){
				for (int timeL = 0; timeL < maxTimeLabel; timeL++){
					double mcTime = timeL*recordStep;
					sf_tDomAveAcc.accum(mcTime, sf_tAcc[domInt][timeL]/4.0);	
					sf_tNDomAveAcc.accum(mcTime, sf_tAcc[nonDomInt][timeL]/4.0);	
					
				}
//			}

			repNo += 1;
			params.set("Reps", repNo);
			writeStDOtoFile();
		}
	}

	void writeStDOtoFile(){
		String message1 = "#Glauber Monte Carlo run: S vs t for several k values. Disorder to Stripes Early times.";
		StringBuffer sb = new StringBuffer();
		String fileName;
				double kRValue = 2*PI*sim.R*2/sim.L;
				sb.append("# k value = ");
				sb.append(kRValue);
				String message2 = sb.toString();
				StringBuffer fileBuf = new StringBuffer();
				
				fileBuf.append(params.sget("Data Dir")); fileBuf.append(File.separator);
				fileBuf.append("d");fileName = fileBuf.toString();
				FileUtil.initFile(fileName, params, message1, message2);
				FileUtil.printAccumToFile(fileName, sf_tDomAveAcc);

				fileBuf.deleteCharAt(fileBuf.length()-1);
				fileBuf.append("n");fileName = fileBuf.toString();
				FileUtil.initFile(fileName, params, message1, message2);
				FileUtil.printAccumToFile(fileName, sf_tNDomAveAcc);
	}
	
	
	private void initialize(){
		sim = new IsingLR(params);
		sim.randomizeField(params.fget("Initial magnetization"));	
		dx = 1;
		clearFile = true;
		sf_tDomAveAcc.clear();
		sf_tNDomAveAcc.clear();
		sf_tDomAveAcc.enableErrorBars(true);
		sf_tNDomAveAcc.enableErrorBars(true);		
	}
	

}

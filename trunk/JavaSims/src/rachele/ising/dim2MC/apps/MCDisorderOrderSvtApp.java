package rachele.ising.dim2MC.apps;

import static java.lang.Math.PI;
import static java.lang.Math.abs;
import static scikit.util.Utilities.format;
import java.awt.Color;
import java.io.File;
//import java.io.File;
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

public class MCDisorderOrderSvtApp extends Simulation{

	Grid grid = new Grid("Long Range Ising Model");
	Grid sfGrid = new Grid("SF");
	Plot sftAvePlot = new Plot("sf_t ave plot");
	Plot sftPlot = new Plot("sf_t plot"); 
	int dx;
	IsingLR sim;
	public FourierTransformer fft;
	double [] sFactor;
	int accNo = 10;
	int [][] sfLabel = new int [accNo][4];
	Accumulator [] sf_tAveAcc = new Accumulator [accNo];
	Accumulator [] sf_tAcc = new Accumulator [accNo];
	boolean clearFile;
	double [] sfTimeArray;

	public static void main(String[] args) {
		new Control(new MCDisorderOrderSvtApp(), "MC Disorder to Order");
	}

	public void load(Control c) {

		//c.frameTogether("grids", grid, sfGrid);		
		c.frameTogether("Data", grid, sftPlot, sftAvePlot);
		c.frame(grid);

		params.add("Data Dir",new DirectoryValue("/Users/erdomi/data/lraim/stripeToClumpInvestigation/mcResults/conservedOP_DO/JR=1/testRuns"));
		params.addm("Dynamics", new ChoiceValue("Kawasaki Glauber", "Ising Glauber", "Kawasaki Metropolis",  "Ising Metropolis"));
		params.add("Random seed", 0);
		params.add("L", 256);
		params.add("R", 23);//1<<6);
		params.addm("Jump Range", 16);
		params.add("Initial magnetization", 0.0);
		params.addm("T", 0.05);
		params.addm("J", -1.0);
		params.addm("h", 0.0);
		params.addm("dt", 0.250);
		params.addm("maxTime", 10.0);
		params.addm("k int", 1);
		params.add("time");
		params.add("magnetization");
		params.add("Lp");
		params.add("Reps");
	}	

	public void animate() {
		sftPlot.setLogScale(false, true);
		sftAvePlot.setLogScale(false,true);
		params.set("time", format(sim.time()));
		params.set("magnetization", format(sim.magnetization()));
		sim.setParameters(params);
		params.set("Lp", sim.L/dx);
		grid.registerData(sim.L/dx, sim.L/dx, sim.getField(dx));
		sfGrid.registerData(sim.L/dx, sim.L/dx, sFactor);
		
		for (int i = 0; i < accNo; i ++){
			StringBuffer sb = new StringBuffer();sb.append("s(t) "); sb.append(i);
			float colorChunk = (float)i/(float)accNo;
			Color col = Color.getHSBColor(colorChunk, 1.0f, 1.0f);
			sftPlot.registerLines(sb.toString(), sf_tAcc[i], col);
			sb.append("Ave");
			sftAvePlot.registerLines(sb.toString(),sf_tAveAcc[i], col);
		}
		

	}

	public void clear() {

	}

	public void run() {
		initialize();
		fft = new FourierTransformer((int)(sim.L/dx));
		sFactor = new double [sim.L/dx*sim.L/dx];
		double step = 0.10;
		double maxTime = params.fget("maxTime");
		int repNo = 0;

		while (true) {
			for (int i = 0; i < accNo; i ++)
				sf_tAcc[i].clear();
			sim.randomizeField(params.fget("Initial magnetization"));
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
					sFactor = fft.calculate2DSF(sim.getField(dx), false, false);
					collect(sFactor);
					recordStep += step;
					recordInt +=1;
				}
			}	
			repNo += 1;
			params.set("Reps", repNo);
			writeToFile();
		}
	}

		public void initialize(){
			sim = new IsingLR(params);
			sim.randomizeField(params.fget("Initial magnetization"));		
			dx = 1;
			clearFile = true;
			int kInt = params.iget("k int");
			for (int i = 0; i < accNo; i++){
				sf_tAcc[i] = new Accumulator();
				sf_tAveAcc[i] = new Accumulator(); sf_tAveAcc[i].enableErrorBars(true);
				sfLabel[i][0] = kInt*i;
				sfLabel[i][1] = kInt*i*sim.L;
				sfLabel[i][2] = (sim.L-kInt*i)%sim.L;
				sfLabel[i][3] = (((sim.L-kInt*i)%sim.L)*sim.L);
				System.out.println("sfLabel " + i + " = " + sfLabel[i][0] + " kRValue hor = " + kRvalue(sfLabel[i][0]));
			}
		}

		private void collect(double [] sFactor){
			for (int i = 0; i < accNo; i ++){
				for(int j=0; j < 4; j++){
					sf_tAveAcc[i].accum(sim.time(),sFactor[sfLabel[i][j]]);			
					sf_tAcc[i].accum(sim.time(),sFactor[sfLabel[i][j]]);
				}
			}
		}
		
		private double kRvalue(int label){
			double kRvalue = 2*PI*label*sim.R/(sim.L); 
			return kRvalue;
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
	
	private void writeToFile(){
		String message1 = "#Glauber Monte Carlo run: S vs t for several values of k. Disorder to order quench.  SF averaged in 4 directions.";
		String fileName = params.sget("Data Dir") + File.separator + "f00";
		StringBuffer fileBuffer = new StringBuffer(); fileBuffer.append(fileName);
		for (int i=0; i < accNo; i ++){
			StringBuffer mb = new StringBuffer();
			mb.append("# k value = ");	mb.append(kRvalue(sfLabel[i][0])); 
			fileBuffer.deleteCharAt(fileBuffer.length()-1);fileBuffer.deleteCharAt(fileBuffer.length()-1);
			if(i < 10) fileBuffer.append("0");
			fileBuffer.append(i); fileName = fileBuffer.toString();
			String message2 = mb.toString();
			FileUtil.initFile(fileName, params, message1, message2);		
			FileUtil.printAccumToFile(fileName, sf_tAveAcc[i]);
		}
		
	}


}


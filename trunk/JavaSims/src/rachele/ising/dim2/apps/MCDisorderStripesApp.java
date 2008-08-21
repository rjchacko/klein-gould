package rachele.ising.dim2.apps;

import static java.lang.Math.PI;
import static scikit.util.Utilities.format;
import java.awt.Color;
import java.io.File;
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

public class MCDisorderStripesApp extends Simulation{
	Grid grid = new Grid("Long Range Ising Model");
	Grid sfGrid = new Grid("SF");
	Plot sfkPlot = new Plot("sf_k plot");
	Plot sftPlot = new Plot("sf_t plot"); 
	int dx;
	IsingLR sim;
	public FourierTransformer fft;
	double [] sFactor;
	Accumulator sf_kAcc; 

	int accNo = 6;
	int [][][] sfLabel = new int [2][4][accNo];//Hor, Vert; four corners to be averaged
	Accumulator [] sf_tDomAveAcc = new Accumulator [accNo]; //Average sf for dominant direction (H or V)
	Accumulator [] sf_tNDomAveAcc = new Accumulator [accNo];//Average sf for non-dominant direction (H or V)
	double [][][] sf_tAcc;// 1st label: 0 = Horizontal direction, 1 = vertical direction; 2nd label -> accNo; 3rd label -> time Label
	Accumulator sf_kTheoryAcc, sf_kTheory2Acc, sf_tTheoryAcc, sf_tTheory2Acc;
    boolean clearFile;
   	double [] sfTimeArray;
	
	public static void main(String[] args) {
		new Control(new MCDisorderStripesApp(), "MC Disorder -> Stripes");
	}
	
	public void load(Control c) {
		
		c.frame(grid);
		c.frame(sftPlot);
		params.add("Data Dir",new DirectoryValue("/home/erdomi/data/lraim/stripeToClumpInvestigation/mcResults/DS/testRuns"));
		params.addm("Dynamics", new ChoiceValue("Ising Glauber","Kawasaki Glauber", "Kawasaki Metropolis",  "Ising Metropolis"));
		params.add("Random seed", 0);
		params.add("L", 1<<7);
		params.add("R", 50);//1<<6);
		params.add("Initial magnetization", 0.0);
		params.addm("T", 0.096548444);
		params.addm("J", -1.0);
		params.addm("h", 0.0);
		params.addm("ky", 2);
		params.addm("dkx", 1);
		params.addm("dt", 1/(double)(1<<3));
		params.addm("maxTime", 10.0);
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
		for (int i = 0; i < accNo; i ++){
			StringBuffer sb1 = new StringBuffer();sb1.append("s(t) Dominant Direction"); sb1.append(i);
			StringBuffer sb2 = new StringBuffer();sb2.append("s(t) NonDominant Direction"); sb2.append(i);
			float colorChunk = (float)i/(float)accNo;
			Color col = Color.getHSBColor(colorChunk, 1.0f, 1.0f);
			sftPlot.registerLines(sb1.toString(), sf_tDomAveAcc[i], col);
			sftPlot.registerLines(sb2.toString(), sf_tNDomAveAcc[i], col);
		}
	}
	
	public void clear() {
	}
	
	public void run() {
		initialize();
		fft = new FourierTransformer((int)(sim.L/dx));
		sFactor = new double [sim.L/dx*sim.L/dx];
		double maxTime = params.fget("maxTime");
		double recordStep = 1/4.0;
		double mod = maxTime%recordStep;
		int maxTimeLabel = (int)(maxTime/recordStep);
		System.out.println("Record Step = " + recordStep + " mod (must be zero) = " + mod);
		sf_tAcc = new double [2][accNo][maxTimeLabel];
		int repNo = 0;

		while (true) {
			//must clear all arrays because we will add 4 contributions to average.
			for (int i = 0; i < 2; i++){
				for (int j = 0; j < accNo; j++){
					for (int k = 0; k < maxTimeLabel; k ++){
						sf_tAcc[i][j][k] = 0.0;						
					}
				}
			}

			sim.randomizeField(params.fget("Initial magnetization"));
			sim.restartClock();
			int timeLabel =0;
			while (timeLabel < maxTimeLabel){
				Job.animate();
				if (sim.time()% recordStep == 0){
					sFactor = fft.calculate2DSF(sim.getField(dx), false, false);
					for (int i = 0; i < accNo; i ++){
						for (int j = 0; j < 4; j++){
							sf_tAcc[0][i][timeLabel] = sFactor[sfLabel[0][j][i]];//horizontal direction
							sf_tAcc[1][i][timeLabel] = sFactor[sfLabel[1][j][i]];//vertical direction
						}
					}
					//double arrayTime = recordStep*timeLabel;
					//System.out.println("check time match: sim time = " + sim.time() + " array time = " + arrayTime);
					timeLabel += 1;
				}
				sim.step();
			}	
			// find hor or vert stripes
			int domInt,nonDomInt;
			if(sFactor[sfLabel[1][0][0]]>sFactor[sfLabel[0][0][0]]){
				domInt = 1; nonDomInt = 0;
				System.out.println("horizontal stripes");
			}else{
				domInt = 0; nonDomInt = 1;
				System.out.println("vertical stripes");
			}
			// record in dom and non-dominant
			for (int i = 0; i < accNo; i++){
				for (int timeL = 0; timeL < maxTimeLabel; timeL++){
					double mcTime = timeL*recordStep;
					sf_tDomAveAcc[i].accum(mcTime, sf_tAcc[domInt][i][timeL]/4.0);	
					sf_tNDomAveAcc[i].accum(mcTime, sf_tAcc[nonDomInt][i][timeL]/4.0);	
				}
			}
			
			repNo += 1;
			params.set("Reps", repNo);
			writeStDOtoFile();
		}
	}

	private void writeStDOtoFile(){
		String message1 = "#Glauber Monte Carlo run: S vs t for several k values. Disorder to Stripes Early times.";
		StringBuffer sb = new StringBuffer();
		String fileName;
			for (int i = 0; i < accNo; i ++){
				double kRValue = 2*PI*sim.R*sfLabel[1][0][i]/sim.L;
				sb.append("# k value = ");
				sb.append(kRValue);
				String message2 = sb.toString();
				StringBuffer fileBuf = new StringBuffer();
				
				fileBuf.append(params.sget("Data Dir")); fileBuf.append(File.separator);
				fileBuf.append("d");fileBuf.append(i);fileName = fileBuf.toString();
				FileUtil.initFile(fileName, params, message1, message2);
				FileUtil.printAccumToFile(fileName, sf_tDomAveAcc[i]);

				fileBuf.deleteCharAt(fileBuf.length()-1);fileBuf.deleteCharAt(fileBuf.length()-1);
				fileBuf.append("n");fileBuf.append(i);fileName = fileBuf.toString();
				FileUtil.initFile(fileName, params, message1, message2);
				FileUtil.printAccumToFile(fileName, sf_tNDomAveAcc[i]);
			}

	}
	
	
	private void initialize(){
		sim = new IsingLR(params);
		sim.randomizeField(params.fget("Initial magnetization"));	
		int ky = params.iget("ky");
		int dkx = params.iget("dkx");
		int index, x, y;
		for (int i =0; i < accNo; i++){
			sf_tDomAveAcc[i] = new Accumulator(); sf_tDomAveAcc[i].enableErrorBars(true);	
			sf_tNDomAveAcc[i] = new Accumulator(); sf_tNDomAveAcc[i].enableErrorBars(true);	
			int kx = i*dkx;
			//Horizontal labels
			int size = sim.L*sim.L;
			y = ky; x = kx;//Up, right
			index = sim.L*y+x; 
			sfLabel[0][0][i] = index;
			x = sim.L - kx;//Up, left	
			index = sim.L*y+x; sfLabel[0][1][i] = index%size;
			y = sim.L-ky; x = kx;//Down, right
			index = sim.L*y+x; sfLabel[0][2][i] = index%size;
			x = sim.L - kx;//Up, left	
			index = sim.L*y+x; sfLabel[0][3][i] = index%size;
			//Vertical labels
			x = ky; y = kx;//Right, Up
			index = sim.L*y+x; sfLabel[1][0][i] = index%size;
			y = sim.L - kx;//Right, Down	
			index = sim.L*y+x; sfLabel[1][1][i] = index%size;
			x = sim.L-ky; y = kx;//Left, Up
			index = sim.L*y+x; sfLabel[1][2][i] = index%size;
			y = sim.L - kx;//Left, Down
			index = sim.L*y+x; sfLabel[1][3][i] = index%size;
		
			if(sfLabel [1][0][i] > size)System.out.println(i + " 0 out");
			if(sfLabel [1][1][i] > size)System.out.println(i + " 1 out");
			if(sfLabel [1][2][i] > size)System.out.println(i + " 2 out");
			if(sfLabel [1][3][i] > size)System.out.println(i + " 3 out");
		}		
	
		dx = 1;
		clearFile = true;
	}

}

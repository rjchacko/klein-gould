package rachele.ising.dim2MC.apps;

import static scikit.util.Utilities.format;
//import java.awt.Color;
//import java.awt.Color;
import java.io.DataInputStream;
import java.io.EOFException;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import rachele.ising.dim2MC.EnergyMetric;
import rachele.ising.dim2MC.IsingLR;
import rachele.util.FileUtil;
import rachele.util.FourierTransformer;
//import scikit.dataset.Accumulator;
import scikit.dataset.Accumulator;
import scikit.graphics.dim2.Grid;
//import scikit.graphics.dim2.Plot;
import scikit.jobs.Control;
import scikit.jobs.Job;
import scikit.jobs.Simulation;
import scikit.jobs.params.ChoiceValue;
import scikit.jobs.params.DirectoryValue;
//import scikit.jobs.params.FileValue;

public class EnergyMetricApp extends Simulation{
	
	IsingLR sim;
	EnergyMetric em;
    FourierTransformer ft;
	
	Grid grid = new Grid("Long Range Ising Model");
//	Plot metricPlot = new Plot("inverseMetric");
//	Plot metric0Plot = new Plot("metric0");
	
	int dx;
	int recordStep = 0;
    boolean clearFile;
    boolean takeAverages;
    String outputFileName, outputFileName0, sfFileNameH, sfFileNameV, sfFileNameDiff;
    
   	Accumulator eMetricInverseAcc = new Accumulator();
   	Accumulator metric0Acc = new Accumulator();
   	
	public static void main(String[] args) {
		new Control(new EnergyMetricApp(), "Ising Model");
	}
	
	public void load(Control c) {
//		c.frameTogether("grids", grid, metricPlot, metric0Plot);	
		c.frame(grid);	
		params.add("Output Directory",new DirectoryValue("/Users/erdomi/data/ergodicityTest/testruns"));
		params.addm("Dynamics", new ChoiceValue("Kawasaki Glauber", "Ising Glauber", "Kawasaki Metropolis",  "Ising Metropolis"));
		params.addm("init", new ChoiceValue( "Random", "Read From File"));
		params.add("Random seed", 0);
		params.add("L", 128);
		params.add("R", 46);
		params.add("JR", 46);
		params.add("Initial magnetization", 0.0);
		params.addm("T", .21);
		params.addm("J", -1.0);
		params.addm("h", 0.0);
		params.addm("dt", .25);
		params.addm("sf int", 2);
		params.add("time");
		params.add("magnetization");
		params.add("Lp");
		params.add("Reps");
		flags.add("Clear");
		flags.add("Metric0");
		flags.add("t=0");
	}	
	
	public void animate() {
		params.set("time", format(sim.time()));
		params.set("magnetization", format(sim.magnetization()));
		sim.setParameters(params);
		params.set("Lp", sim.L/dx);
		
		grid.registerData(sim.L/dx, sim.L/dx, sim.getField(dx));
//		metricPlot.setAutoScale(true);
//		metric0Plot.setAutoScale(true);
//		metricPlot.registerPoints("eMet", eMetricInverseAcc, Color.blue);
//		metric0Plot.registerPoints("metric0", metric0Acc, Color.red);
		
		if(flags.contains("Clear")){
			eMetricInverseAcc.clear();
			metric0Acc.clear();
		}
		if(flags.contains("t=0")){
			sim.resetTime();
			recordStep = 0;
			em.clearSums();
			FileUtil.deleteFile(outputFileName);
			FileUtil.deleteFile(outputFileName0);
//			FileUtil.deleteFile(sfFileNameH);
//			FileUtil.deleteFile(sfFileNameV);
//			FileUtil.deleteFile(sfFileNameDiff);
			FileUtil.initFile(outputFileName, params);
			FileUtil.initFile(outputFileName0, params);
			
		}
		if(flags.contains("Metric0")){
			double metric0 = em.findMetric0();
			System.out.println("Metric0 = " + metric0);
			metric0Acc.accum(sim.time(), metric0);
		}
		flags.clear();

	}
	
	public void clear() {
		eMetricInverseAcc.clear();
		metric0Acc.clear();
	}
	
	public void run() {
		
		initialize();
		recordStep = 0;
		
		double [] sf = new double [sim.L*sim.L];
		em = new EnergyMetric(sim, params);
		
		double step = 0.0;
		int sfInt = params.iget("sf int");
		
		while (true) {
			sim.step();
			em.calculateMetric();
			if (sim.time() > recordStep){
				double pt = 1.0/em.eMetric;
				double met0 = em.findMetric0();
//				eMetricInverseAcc.accum(sim.time(), pt);
//				metric0Acc.accum(sim.time(), met0);
				FileUtil.printlnToFile(outputFileName, sim.time(), pt);
				FileUtil.printlnToFile(outputFileName0, sim.time(), met0);
				sf = ft.calculate2DSF(sim.getField(1), false, false);
				double hpoint = 0.5*(sf[sfInt] +sf[(sim.L-sfInt)]);
				double vpoint = 0.5*(sf[sfInt*sim.L] +sf[(sim.L-sfInt)*sim.L]);
				double diff = (hpoint - vpoint);
				FileUtil.printlnToFile(sfFileNameH, sim.time(), hpoint);
				FileUtil.printlnToFile(sfFileNameV, sim.time(), vpoint);
				FileUtil.printlnToFile(sfFileNameDiff, sim.time(), diff);
				recordStep += step;
			}
			Job.animate();
		}
	}

	public void writeToFile(String fileName){
		FileUtil.initFile(fileName, params);
		
	}
	
	public void initialize(){
		sim = new IsingLR(params);
		//sim.setField(params.fget("Initial magnetization"));
		sim.randomizeField(params.fget("Initial magnetization"));		
		dx = 1;
		if(params.sget("init") == "Read From File") readInitialConfiguration();
		clearFile = true;
		outputFileName = params.sget("Output Directory") + File.separator + "T" + "test" + ".txt";
		outputFileName0 = params.sget("Output Directory") + File.separator + "0T" + "test" + ".txt";
		sfFileNameH = params.sget("Output Directory") + File.separator + "T" + "test" + "h.txt";
		sfFileNameV = params.sget("Output Directory") + File.separator + "T" + "test" + "v.txt";
		sfFileNameDiff = params.sget("Output Directory") + File.separator + "T" + "test" + "d.txt";
//
//		outputFileName = params.sget("Output Directory") + File.separator + "T" + sim.T + ".txt";
//		outputFileName0 = params.sget("Output Directory") + File.separator + "0T" + sim.T + ".txt";
//		sfFileNameH = params.sget("Output Directory") + File.separator + "T" + sim.T + "h.txt";
//		sfFileNameV = params.sget("Output Directory") + File.separator + "T" + sim.T + "v.txt";
//		sfFileNameDiff = params.sget("Output Directory") + File.separator + "T" + sim.T + "d.txt";
		
		FileUtil.initFile(outputFileName, params);
		FileUtil.initFile(outputFileName0, params);		
		FileUtil.initFile(sfFileNameH, params, "#Horizontal structure factor");			
		FileUtil.initFile(sfFileNameV, params, "#Vertical tructure factor");	
		FileUtil.initFile(sfFileNameDiff, params, "#Abs value of difference of hor & vert");

		ft = new FourierTransformer(sim.L);
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

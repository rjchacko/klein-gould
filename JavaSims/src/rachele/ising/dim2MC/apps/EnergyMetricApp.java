package rachele.ising.dim2MC.apps;

import static scikit.util.Utilities.format;
//import java.awt.Color;
import java.io.DataInputStream;
import java.io.EOFException;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import rachele.ising.dim2MC.EnergyMetric;
import rachele.ising.dim2MC.IsingLR;
import rachele.util.FileUtil;
//import scikit.dataset.Accumulator;
import scikit.graphics.dim2.Grid;
import scikit.graphics.dim2.Plot;
import scikit.jobs.Control;
import scikit.jobs.Job;
import scikit.jobs.Simulation;
import scikit.jobs.params.ChoiceValue;
import scikit.jobs.params.DirectoryValue;
//import scikit.jobs.params.FileValue;

public class EnergyMetricApp extends Simulation{
	
	IsingLR sim;
	EnergyMetric em;
	
	Grid grid = new Grid("Long Range Ising Model");
	Plot metricPlot = new Plot("inverseMetric");
	
	int dx;
    boolean clearFile;
    boolean takeAverages;
    String outputFileName;
	
//   	Accumulator eMetricInverseAcc = new Accumulator();
   	
	public static void main(String[] args) {
		new Control(new EnergyMetricApp(), "Ising Model");
	}
	
	public void load(Control c) {
		c.frameTogether("grids", grid, metricPlot);		
		params.add("Output Directory",new DirectoryValue("/Users/erdomi/data/ergodicityTest/R64/"));
		params.addm("Dynamics", new ChoiceValue("Ising Glauber","Kawasaki Glauber", "Kawasaki Metropolis",  "Ising Metropolis"));
		params.addm("init", new ChoiceValue( "Random", "Read From File"));
		params.add("Random seed", 0);
		params.add("L", 1<<9);
		params.add("R", 1<<6);
		params.add("Initial magnetization", 0.0);
		params.addm("T", 1.00001);
		params.addm("J", 1.0);
		params.addm("h", 0.0);
		params.addm("dt", .25);
		params.add("time");
		params.add("magnetization");
		params.add("Lp");
		params.add("Reps");
		flags.add("Clear");
	}	
	
	public void animate() {
		params.set("time", format(sim.time()));
		params.set("magnetization", format(sim.magnetization()));
		sim.setParameters(params);
		params.set("Lp", sim.L/dx);
		
		grid.registerData(sim.L/dx, sim.L/dx, sim.getField(dx));
		metricPlot.setAutoScale(true);
//		metricPlot.registerPoints("eMet", eMetricInverseAcc, Color.blue);
		
		if(flags.contains("Clear")){
//			eMetricInverseAcc.clear();
			flags.clear();
		}

	}
	
	public void clear() {
//		eMetricInverseAcc.clear();
	}
	
	public void run() {
		initialize();
		em = new EnergyMetric(sim, params);
		double step = 0.10;
		int recordStep = 0;
		while (true) {
			sim.step();
			em.calculateMetric();
			if (sim.time() > recordStep){
//				eMetricInverseAcc.accum(sim.time(), 1.0/em.eMetric);
				double pt = 1.0/em.eMetric;
				FileUtil.printlnToFile(outputFileName, sim.time(), pt);
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
		outputFileName = params.sget("Output Directory") + "eT" + sim.T + ".txt";
		FileUtil.initFile(outputFileName, params);
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
	
	public void recordSfDataToFile(double [] data){
		String file0 = "../../../research/javaData/stripeToClumpInvestigation/monteCarloData/s0";
		String file1 = "../../../research/javaData/stripeToClumpInvestigation/monteCarloData/s1";
		String file2 = "../../../research/javaData/stripeToClumpInvestigation/monteCarloData/s2";
		String file3 = "../../../research/javaData/stripeToClumpInvestigation/monteCarloData/s3";
		String file4 = "../../../research/javaData/stripeToClumpInvestigation/monteCarloData/s4";
		String file5 = "../../../research/javaData/stripeToClumpInvestigation/monteCarloData/s5";
		if (clearFile){
			FileUtil.deleteFile(file0);
			FileUtil.deleteFile(file1);
			FileUtil.deleteFile(file2);
			FileUtil.deleteFile(file3);
			FileUtil.deleteFile(file4);
			FileUtil.deleteFile(file5);
			clearFile = false;
		}
		int ky = 2;
		int Lp = sim.L/dx;
		FileUtil.printlnToFile(file0, sim.time(), data[ky*Lp]*data[ky*Lp]);
		FileUtil.printlnToFile(file1, sim.time(), data[ky*Lp+1]*data[ky*Lp+1]);
		FileUtil.printlnToFile(file2, sim.time(), data[ky*Lp+2]*data[ky*Lp+2]);
		FileUtil.printlnToFile(file3, sim.time(), data[ky*Lp+3]*data[ky*Lp+3]);
		FileUtil.printlnToFile(file4, sim.time(), data[ky*Lp+4]*data[ky*Lp+4]);
		FileUtil.printlnToFile(file5, sim.time(), data[ky*Lp+5]*data[ky*Lp+5]);
		//System.out.println("data written for time = " + ising.time());
	}	
	
}

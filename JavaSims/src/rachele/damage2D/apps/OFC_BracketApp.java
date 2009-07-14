package rachele.damage2D.apps;

import java.io.File;

import rachele.damage2D.OFC_Lattice;
import rachele.util.FileUtil;
import scikit.dataset.Accumulator;
import scikit.dataset.Histogram;
import scikit.graphics.dim2.Grid;
import scikit.jobs.Control;
import scikit.jobs.Job;
import scikit.jobs.Simulation;
import scikit.jobs.params.DirectoryValue;
import scikit.util.Utilities;

public class OFC_BracketApp extends Simulation{

	int cg_dt;

	Grid grid = new Grid("Lattice");
	Grid cgGrid = new Grid(" CG grid");
	OFC_Lattice ofc;
	Accumulator cgMetricAcc = new Accumulator();
	Accumulator sizeStore = new Accumulator();  //To store data at each plate update within a bracket
	Accumulator stressMetricStore = new Accumulator();  //To store data at each plate update within a bracket
	Histogram sizeHist = new Histogram(1);
	String iMetBracketFile;
	String sizeBracketFile;


	public static void main(String[] args) {
		new Control(new OFC_BracketApp(), "OFC Model Bracket App");
	}

	public void load(Control c) {
		c.frameTogether("Grids", grid, cgGrid);
		params.add("Data Dir",new DirectoryValue("/Users/erdomi/data/damage/testRuns"));
		params.addm("Random Seed", 1);
		params.addm("CG size", 32);
		params.addm("dx", 8);
		params.addm("Coarse Grained dt", 500);
		params.addm("Equilibration Updates", 1000000);
		params.addm("Bracket Min", 600000);
		params.addm("Bracket Max", 700000);
		params.addm("R", 16);// 0 -> fully connected
		params.addm("Residual Stress", 0.625);
		params.addm("Dissipation Param", 0.05);
		params.addm("Res. Max Noise", 0.125);
		params.addm("Lower Cutoff", 1);
		params.add("L");
		params.add("Time");
		params.add("Plate Updates");
	}

	public void animate() {
		grid.registerData(ofc.L, ofc.L, ofc.stress);
		cgGrid.registerData(ofc.Lp, ofc.Lp, ofc.epicenterCount);
		params.set("Time", Utilities.format(ofc.cg_time));
		params.set("Plate Updates", ofc.plateUpdates);
	}

	public void clear() {
	}

	public void run() {
		ofc = new OFC_Lattice(params);
		initFiles();
		cg_dt = params.iget("Coarse Grained dt");
		int lowerBracket = params.iget("Bracket Min");
		int upperBracket = params.iget("Bracket Max");

		//equilibrate
		ofc.initEquilibrate(params.iget("Equilibration Updates"));
		while (ofc.plateUpdates < 0){
			ofc.equilibrate();
			Job.animate();
		}

		while(true){

			ofc.step();
			int size = ofc.avSize;
			double iMet = ofc.calcInverseMetric();
			if(ofc.cg_time > lowerBracket){
				double reducedTime = ofc.cg_time/cg_dt;
				FileUtil.printlnToFile(iMetBracketFile, ofc.cg_time, iMet, reducedTime, size);
			}

			if(ofc.cg_time > upperBracket) Job.signalStop();

			Job.animate();
		}
	}

	void initFiles(){
		iMetBracketFile = params.sget("Data Dir") + File.separator + "imb.txt";  // to record iverse metric data
		FileUtil.initFile(iMetBracketFile, params, " time (plate updates), stress inverse metric, time/coarse grained time, size of avalanche");	
	}
}

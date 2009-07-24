package rachele.damage2D.apps;

import java.awt.Color;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;

import rachele.damage2D.OFC_Lattice;
import rachele.util.FileUtil;
import scikit.dataset.Accumulator;
import scikit.dataset.DatasetBuffer;
import scikit.dataset.Histogram;
import scikit.graphics.dim2.Geom2D;
import scikit.graphics.dim2.Grid;
import scikit.jobs.Control;
import scikit.jobs.Job;
import scikit.jobs.Simulation;
import scikit.jobs.params.DirectoryValue;
import scikit.util.Utilities;

public class OFC_DelayWriteApp extends Simulation{

	int cg_dt;
	Grid grid = new Grid("Lattice");
	Grid cgGrid = new Grid(" CG grid");
	Grid cgGridTimeAverage = new Grid("Time ave CG grid");
	Grid plateUpdateGrid = new Grid("Plate Update grid");
	OFC_Lattice ofc;
	Histogram sizeHist = new Histogram(1);
	String iMetFile;
	String sizeHistFile;
	int maxSize;
	
	Accumulator invStressMetTempAcc = new Accumulator();
	Accumulator cgInvStressMetTempAcc = new Accumulator();
	Accumulator cgInvActivityMetTempAcc = new Accumulator();
	Accumulator maxSizeTempAcc = new Accumulator();
	
	public static void main(String[] args) {
		new Control(new OFC_DelayWriteApp(), "OFC Model");
	}
	
	public void load(Control c) {
		c.frameTogether("Grids", grid, cgGrid, plateUpdateGrid, cgGridTimeAverage);
		params.add("Data Dir",new DirectoryValue("/Users/erdomi/data/damage/testRuns"));
		params.addm("Random Seed", 1);
		params.addm("CG size", 16);
		params.addm("dx", 1);
		params.addm("Coarse Grained dt (PU)", 100);
		params.addm("Equilibration Updates", 500000);
		params.addm("Max CG Time", 1000000);
		params.addm("Data points per write", 1000);
		params.addm("R", 4);// 0 -> fully connected
		params.addm("Residual Stress", 0.625);
		params.addm("Dissipation Param", 0.1);
		params.addm("Res. Max Noise", 0.125);
		params.addm("Lower Cutoff", 1);
		params.add("L");
		params.add("CG Time");
		params.add("Size");
		params.add("Plate Updates");
	}
	
	public void animate() {
		grid.registerData(ofc.L, ofc.L, ofc.stress);
		cgGrid.registerData(ofc.Lp, ofc.Lp, ofc.epicenterCount);
		cgGridTimeAverage.registerData(ofc.Lp, ofc.Lp, ofc.CG_ActivityTimeAve);
		plateUpdateGrid.registerData(ofc.L, ofc.L, ofc.plateUpdateFailLocations);
		
		grid.clearDrawables();
		double radius = 1.0/(2.0*ofc.L);
		double failSite_y = ((double)(ofc.epicenterSite/ofc.L))/ofc.L + radius;
		double failSite_x = ((double)(ofc.epicenterSite%ofc.L))/ofc.L + radius;
		grid.addDrawable(
				Geom2D.circle(failSite_x, failSite_y, radius, Color.GREEN));
		params.set("CG Time", Utilities.format(ofc.cg_time));
		params.set("Plate Updates", ofc.plateUpdates);
		params.set("Size", ofc.avSize);
	}

	public void clear() {
	}

	public void run() {
		ofc = new OFC_Lattice(params);
		initFiles();
		cg_dt = params.iget("Coarse Grained dt (PU)");
		double nextAccumTime = 0;
		int maxTime = params.iget("Max CG Time");
		int dataPointsPerWrite = params.iget("Data points per write");
		
		maxSize = 0;
		int dataPointsCount = 0;
		//equilibrate
		ofc.initEquilibrate(params.iget("Equilibration Updates"));
		while (ofc.plateUpdates < 0){
			ofc.equilibrate();
			Job.animate();
		}
		
		while(true){
			ofc.step();

			int size = ofc.avSize;
			sizeHist.accum(size);
			if (size > maxSize) {
				maxSize = size;
			}
			double stressMetric = ofc.calcStressMetric();

			if(ofc.plateUpdates > nextAccumTime){ //Accumulate data to be written
				Job.animate();
				//stress metric
				double inverseStressMetric = 1.0/stressMetric;
				invStressMetTempAcc.accum(ofc.plateUpdates, inverseStressMetric);
				//CG stress metric
				double cgInverseStressMetric = 1.0/ofc.calcCG_stressMetric();				
				cgInvStressMetTempAcc.accum(ofc.cg_time, cgInverseStressMetric);
				//CG activity metric
				double cgInverseActivityMetric = 1.0/ofc.calcCG_activityMetric();
				cgInvActivityMetTempAcc.accum(ofc.cg_time, cgInverseActivityMetric);
				//maxSize
				maxSizeTempAcc.accum(ofc.cg_time, maxSize);
				
				nextAccumTime += cg_dt;
				dataPointsCount += 1;
				maxSize = 0;

			}
			
//			System.out.println("no data pts = " + dataPointsCount + " cg time = " + ofc.cg_time + " next accum time = " + nextAccumTime);
			
			if(dataPointsCount > dataPointsPerWrite -1){ //write accumulated data
				FileUtil.initFile(sizeHistFile, params, "avalanch size histogram");
				FileUtil.printHistToFile(sizeHistFile, sizeHist);
				printAccumsToFile(iMetFile, invStressMetTempAcc, cgInvStressMetTempAcc, cgInvActivityMetTempAcc, maxSizeTempAcc);
				//				FileUtil.printlnToFile(iMetFile, ofc.plateUpdates, inverseStressMetric, ofc.cg_time, cgInverseActivityMetric, cgInverseStressMetric, maxSize);
			
				//restart counts
				maxSize = 0;
				dataPointsCount = 0;
				
				//clear temp accumulators
				invStressMetTempAcc.clear();
				cgInvStressMetTempAcc.clear();
				cgInvActivityMetTempAcc.clear();
				maxSizeTempAcc.clear();

			}

		
			if(ofc.plateUpdates > maxTime) Job.signalStop();
				

		}
	}
	
	static void printAccumsToFile(String fileName, Accumulator acc1, Accumulator acc2, Accumulator acc3, Accumulator acc4){
		DatasetBuffer data1 = acc1.copyData();
		DatasetBuffer data2 = acc2.copyData();
		DatasetBuffer data3 = acc3.copyData();
		DatasetBuffer data4 = acc4.copyData();
		int size = data1.size();
		try{
			File file = new File(fileName);
			PrintWriter pw = new PrintWriter(new FileWriter(file, true), true);
			for (int i = 0; i < size; i++){
				pw.println(data1.x(i) + " " + data1.y(i) + " " + data2.x(i) + " " + data2.y(i) + " " + data3.y(i) + " " + data4.y(i)) ;
			}
		} catch (IOException ex){
			ex.printStackTrace();
		}
	}
	
	void initFiles(){
		iMetFile = params.sget("Data Dir") + File.separator + "im.txt";  // to record iverse metric data
		FileUtil.initFile(iMetFile, params, " time (plate updates), stress inverse metric, time/coarse grained time, coarse grained activity metric, coarse grained stress metric, size of avalanche");
//		sizeFile = params.sget("Data Dir") + File.separator + "s.txt";   //to record size vs time data
//		FileUtil.initFile(sizeFile, params, "avalanch size vs time");
		sizeHistFile = params.sget("Data Dir") + File.separator + "sh.txt";	//to record size histogram data
//		cgInvMetricFile = params.sget("Data Dir") + File.separator + "cg.txt";
//		FileUtil.initFile(cgInvMetricFile, params, "Coarse Grained Inverse Metric");
	}

}

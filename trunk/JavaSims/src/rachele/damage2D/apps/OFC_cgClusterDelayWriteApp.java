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

public class OFC_cgClusterDelayWriteApp extends Simulation{

	int dt;
	Grid grid = new Grid("Lattice");
	Grid cgGrid = new Grid(" CG grid");
	Grid cgGridTimeAverage = new Grid("Time ave CG grid");
	Grid plateUpdateGrid = new Grid("Plate Update grid");
	Grid cgGridActTimeAverage = new Grid("Time ave CG Act Grid");
	Grid cgGridTimeAveDev = new Grid("Time ave CG Dev grid");
	Grid cgGridActTimeAveDev = new Grid("Time ave CG Act Dev Grid");
	OFC_Lattice ofc;
	Histogram sizeHist = new Histogram(1);
	String iMetFile;
	String sizeHistFile;
	String sizeFile;
	String specialBoxFile;
	String overallBoxFile;
	int maxSize;
	int specialBox;
	int lowerCgBound;
	double [] cgCounter;
	
	Accumulator invStressMetTempAcc = new Accumulator();
	Accumulator cgInvStressMetTempAcc = new Accumulator();
	Accumulator cgInvActivityMetTempAcc = new Accumulator();
	Accumulator maxSizeTempAcc = new Accumulator();
	Accumulator invSizeActMetTempAcc = new Accumulator();
	Accumulator	cgMaxSizeActBoxTempAcc = new Accumulator();
	Accumulator	cgMaxSizeActLocationTempAcc = new Accumulator();
	Accumulator timeAveForMaxCGSizeActTempAcc = new Accumulator();
	Accumulator specialBoxTempAcc = new Accumulator();
	Accumulator specialBoxActTempAcc = new Accumulator();
	Accumulator specialBoxSizeActTempAcc = new Accumulator();
	Accumulator specialBoxCGStressTempAcc = new Accumulator();
	Accumulator overallCgCountTempAcc = new Accumulator();
	
	
	
	public static void main(String[] args) {
		new Control(new OFC_cgClusterDelayWriteApp(), "OFC Model- clusters");
	}
	
	public void load(Control c) {
		c.frameTogether("Grids", grid, cgGrid, plateUpdateGrid, cgGridTimeAveDev, cgGridActTimeAveDev);
		params.add("Data Dir",new DirectoryValue("/Users/erdomi/data/damage/testRuns"));
		params.addm("Random Seed", 1);
		params.addm("CG size", 8);
		params.addm("dx", 16);
		params.addm("Coarse Grained dt (PU)", 1);
		params.addm("Equilibration Updates", 500000);
		params.addm("Max PU", 1000000);
		params.addm("Data points per write", 100);
		params.addm("R", 8);// 0 -> fully connected
		params.addm("Residual Stress", 0.625);
		params.addm("Dissipation Param", 0.05);
		params.addm("Res. Max Noise", 0.125);
		params.addm("Lower Cutoff", 1);
		params.addm("Lower CG Bound fraction", 0.07);
		params.addm("Special Box", 0);
		params.addm("Next Pause", 400000);
		params.add("L");
		params.add("CG Time");
		params.add("Size");
		params.add("Plate Updates");
	}
	
	public void animate() {
		grid.registerData(ofc.L, ofc.L, ofc.stress);
		cgGrid.registerData(ofc.Lp, ofc.Lp, ofc.epicenterCount);
		cgGridTimeAveDev.registerData(ofc.Lp, ofc.Lp, ofc.CG_SizeActTimeAveDev);
		cgGridActTimeAveDev.registerData(ofc.Lp, ofc.Lp, ofc.CG_ActivityTimeAveDev);
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
		dt = params.iget("Coarse Grained dt (PU)");
		double nextAccumTime = 0;
		int maxTime = params.iget("Max PU");
		int dataPointsPerWrite = params.iget("Data points per write");
		specialBox = params.iget("Special Box");
		lowerCgBound = (int)((double)(ofc.N)*params.fget("Lower CG Bound fraction"));
		System.out.println("Lower bound = " + lowerCgBound);
		FileUtil.printlnToFile(overallBoxFile, "# lower cg bound = " + lowerCgBound);
		clearTempAccs();
		
		maxSize = 0;
		int dataPointsCount = 0;
		//equilibrate
		ofc.initEquilibrate(params.iget("Equilibration Updates"));
		while (ofc.plateUpdates < 0){
			ofc.equilibrate();
			Job.animate();
		}
		
		cgCounter = new double [ofc.Np];
		for (int i = 0; i < ofc.Np; i++) cgCounter [i] = -1;
		
		while(true){

			ofc.step();

			int size = ofc.avSize;
			sizeHist.accum(size);
			if (size > maxSize) {
				maxSize = size;
			}
			double stressMetric = ofc.calcStressMetricPU();
//			double stressMetric = ofc.calcStressMetric();

			if(ofc.plateUpdates > nextAccumTime){ //Accumulate data to be written
				Job.animate();
				
				int overallCgCount = 0;
				
				//maxSize
				maxSizeTempAcc.accum(ofc.plateUpdates, maxSize);
		
				int maxEpicenter=0;
				int loc = 0;
				for(int i = 0; i < ofc.Np; i++){
					if (ofc.epicenterSize[i] > maxEpicenter){
						maxEpicenter = ofc.epicenterSize[i];
						loc = i;
					}
					if (ofc.epicenterSize[i] > lowerCgBound){
						System.out.println("Box " + i + " has size " + ofc.epicenterSize[i]);
						cgCounter [i] += 1;
						if(cgCounter[i] > 0){
							overallCgCount += 1;
						}
					}
					
				}
				
				cgMaxSizeActBoxTempAcc.accum(ofc.cg_time, maxEpicenter);
				cgMaxSizeActLocationTempAcc.accum(ofc.cg_time, loc);
				specialBoxTempAcc.accum((double)ofc.plateUpdates, (double)ofc.epicenterSize[specialBox]);
//				if(loc==specialBox) System.out.println("special box size = " + ofc.epicenterSize[loc] + " " +loc + " mec " + maxEpicenter + " at " + ofc.plateUpdates);
				specialBoxActTempAcc.accum(ofc.plateUpdates, (double)ofc.epicenterCount[specialBox]);
				specialBoxSizeActTempAcc.accum(ofc.plateUpdates, (double)ofc.epicenterSize[specialBox]);
				specialBoxCGStressTempAcc.accum(ofc.plateUpdates, (double)ofc.CG_Stress[specialBox]);
				overallCgCountTempAcc.accum(ofc.plateUpdates, overallCgCount);
				
				timeAveForMaxCGSizeActTempAcc.accum(ofc.cg_time, ofc.CG_SizeActTimeAve[loc]);
				if(dataPointsCount ==0) System.out.println("max loc = " + maxEpicenter + " " + ofc.CG_SizeActTimeAve[loc]);
				
		
				
				//stress metric
				double inverseStressMetric = 1.0/stressMetric;
				invStressMetTempAcc.accum(ofc.plateUpdates, inverseStressMetric);
				//CG stress metric
				double cgInverseStressMetric = 1.0/ofc.calcCG_stressMetricPU();				
//				double cgInverseStressMetric = 1.0/ofc.calcCG_stressMetric();				
				cgInvStressMetTempAcc.accum(ofc.cg_time, cgInverseStressMetric);
				//CG activity metric
				double cgInverseActivityMetric = 1.0/ofc.calcCG_activityMetricPU();
//				double cgInverseActivityMetric = 1.0/ofc.calcCG_activityMetric();
				cgInvActivityMetTempAcc.accum(ofc.cg_time, cgInverseActivityMetric);
				// size activity metric
				double cgInverseSizeActMet = 1.0/ofc.calcCG_sizeActMetricPU();
//				double cgInverseSizeActMet = 1.0/ofc.calcCG_sizeActMetric();
				invSizeActMetTempAcc.accum(ofc.cg_time, cgInverseSizeActMet);
				

				
				nextAccumTime += dt;
				dataPointsCount += 1;
				maxSize = 0;

			}
			
//			System.out.println("no data pts = " + dataPointsCount + " cg time = " + ofc.cg_time + " next accum time = " + nextAccumTime);
			
			if(dataPointsCount > dataPointsPerWrite -1){ //write accumulated data
				FileUtil.initFile(sizeHistFile, params, "avalanch size histogram");
				FileUtil.printHistToFile(sizeHistFile, sizeHist);
				printAccumsToFile(iMetFile, invStressMetTempAcc, cgInvStressMetTempAcc, cgInvActivityMetTempAcc, invSizeActMetTempAcc);
				printAccumsToFile(sizeFile, maxSizeTempAcc, cgMaxSizeActBoxTempAcc, timeAveForMaxCGSizeActTempAcc, cgMaxSizeActLocationTempAcc);
				printAccumsToFile(specialBoxFile, specialBoxTempAcc, specialBoxCGStressTempAcc, specialBoxActTempAcc, specialBoxSizeActTempAcc);
				FileUtil.printAccumToFileNoErrorBars(overallBoxFile, overallCgCountTempAcc);
				//restart counts
				maxSize = 0;
				dataPointsCount = 0;
				
				//clear temp accumulators
				clearTempAccs();


			}

		
			int nextPause = params.iget("Next Pause");
			if(ofc.plateUpdates == nextPause) Job.signalStop();
			
			if(ofc.plateUpdates > maxTime) Job.signalStop();
//			if(ofc.plateUpdates == 500000) Job.signalStop();
//			if(ofc.plateUpdates == 500001) Job.signalStop();
				

		}
	}
	
	void clearTempAccs(){
		invStressMetTempAcc.clear();
		cgInvStressMetTempAcc.clear();
		cgInvActivityMetTempAcc.clear();
		maxSizeTempAcc.clear();
		invSizeActMetTempAcc.clear();
		cgMaxSizeActBoxTempAcc.clear();
		cgMaxSizeActLocationTempAcc.clear();
		timeAveForMaxCGSizeActTempAcc.clear();
		specialBoxTempAcc.clear();
		specialBoxActTempAcc.clear();
		specialBoxSizeActTempAcc.clear();
		specialBoxCGStressTempAcc.clear();
		overallCgCountTempAcc.clear();
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
	
	static void printAccumsToFile2(String fileName, Accumulator acc1, Accumulator acc2, Accumulator acc3){
		DatasetBuffer data1 = acc1.copyData();
		DatasetBuffer data2 = acc2.copyData();
		DatasetBuffer data3 = acc3.copyData();
		int size = data1.size();
		try{
			File file = new File(fileName);
			PrintWriter pw = new PrintWriter(new FileWriter(file, true), true);
			for (int i = 0; i < size; i++){
				pw.println(data1.x(i) + " " + data1.y(i) + " " + data2.y(i) + " " + data3.y(i) ) ;
			}
		} catch (IOException ex){
			ex.printStackTrace();
		}
	}
	
	void initFiles(){
		iMetFile = params.sget("Data Dir") + File.separator + "im.txt";  // to record iverse metric data
		FileUtil.initFile(iMetFile, params, " time (plate updates), stress inverse metric, time/coarse grained time, coarse grained activity inverse metric, coarse grained stress inverse metric, cg size act inverse metric");
		sizeFile = params.sget("Data Dir") + File.separator + "s.txt";
		FileUtil.initFile(sizeFile, params, "time (plate updates, fixed),  max avalanche size, cgTime, max size act, time av of max size act, loc of max size act");
		sizeHistFile = params.sget("Data Dir") + File.separator + "sh.txt";	//to record size histogram data
		specialBoxFile = params.sget("Data Dir") + File.separator + "b.txt";
		FileUtil.initFile(specialBoxFile, params, "time (Plate Updates), size-activity");
		overallBoxFile = params.sget("Data Dir") + File.separator + "ab.txt";
		FileUtil.initFile(overallBoxFile, params, "overall box file: counts number of times that a cg box has a reoccurance of SA > lower cg bound", "time (Plate Updates), no of reoccurances");	
	}

}

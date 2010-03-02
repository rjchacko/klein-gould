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
import scikit.graphics.dim2.Plot;
import scikit.jobs.Control;
import scikit.jobs.Job;
import scikit.jobs.Simulation;
import scikit.jobs.params.DirectoryValue;
import scikit.util.Utilities;

public class OFC_DelayWriteApp extends Simulation{

	int dt;
	Grid grid = new Grid("Lattice");
	Grid cgGrid = new Grid(" CG grid");
	Grid cgGridTimeAverage = new Grid("Time ave CG grid");
	Grid plateUpdateGrid = new Grid("Plate Update grid");
	Plot histPlot = new Plot("Stress histogram");
	Plot saHistPlot = new Plot("SA histogram");
	Plot amMetBreakdown = new Plot("Met Breakdown");
	Plot samMetBreakdown = new Plot("SA Met Breakdown");
	OFC_Lattice ofc;
	Histogram sizeHist = new Histogram(1);
	String iMetFile;
	String sizeHistFile;
	String sizeFile;

	int maxSize;
	
	Accumulator invStressMetTempAcc = new Accumulator();
	Accumulator cgInvStressMetTempAcc = new Accumulator();
	Accumulator cgInvActivityMetTempAcc = new Accumulator();
	Accumulator maxSizeTempAcc = new Accumulator();
	Accumulator invSizeActMetTempAcc = new Accumulator();
	Accumulator	cgMaxSizeActBoxTempAcc = new Accumulator();
	Accumulator	cgMaxSizeActLocationTempAcc = new Accumulator();
	Accumulator timeAveForMaxCGSizeActTempAcc = new Accumulator();
	Histogram stressHist = new Histogram(0.001);
	Histogram saHist = new Histogram(0.0001);
	Accumulator amSqAveAccum = new Accumulator ();
	Accumulator amAveSqAccum = new Accumulator ();
	Accumulator samSqAveAccum = new Accumulator ();
	Accumulator samAveSqAccum = new Accumulator ();
	
	public static void main(String[] args) {
		new Control(new OFC_DelayWriteApp(), "OFC Model");
	}
	
	public void load(Control c) {
		c.frameTogether("Grids", grid, cgGrid, plateUpdateGrid, cgGridTimeAverage);
//		c.frame(histPlot);
//		c.frame(saHistPlot);
//		c.frame(amMetBreakdown);
//		c.frame(samMetBreakdown);
		params.add("Data Dir",new DirectoryValue("/Users/erdomi/data/damage/earthquakeOnly/smallTimeScaleTeste100/"));
		params.addm("Random Seed", 1);
		params.addm("CG size", 8);
		params.addm("dx", 2);
		params.addm("Coarse Grained dt (PU)", 1000);
		params.addm("Equilibration Updates", 50000);
		params.addm("Max PU", 1000000);
		params.addm("Data points per write", 1);
		params.addm("R", 16);// 0 -> fully connected
		params.addm("Residual Stress", 0.625);
		params.addm("Dissipation Param", 0.3);
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
		
		
		//Distributions
//		stressHist.clear();
//		saHist.clear();
//		for(int i = 0; i<ofc.N; i++){
//			stressHist.accum(ofc.stressTimeAve[i]);
//		}
//		for(int i = 0; i<ofc.Np; i++){
//			saHist.accum(ofc.CG_SizeActTimeAve[i]);
//		}
//		saHistPlot.registerBars("SA", saHist, Color.BLUE);
//		histPlot.registerBars("Stress", stressHist, Color.RED);
		
		//Metric calc breakdowns
//		amMetBreakdown.registerLines("am sq ave", amSqAveAccum, Color.BLACK);
//		amMetBreakdown.registerLines("am ave sq", amAveSqAccum, Color.BLUE);
//		
//		samMetBreakdown.registerLines("sam sq ave", samSqAveAccum, Color.RED);
//		samMetBreakdown.registerLines("sam ave sq", samAveSqAccum, Color.GREEN);
		
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
			double stressMetric = ofc.calcStressMetricPU();
//			double stressMetric = ofc.calcStressMetric();

			if(ofc.plateUpdates >= nextAccumTime){ //Accumulate data to be written
//				System.out.println("nextAccumTime = " + nextAccumTime);
				Job.animate();
				
				//maxSize
				maxSizeTempAcc.accum(ofc.plateUpdates, maxSize);
				
				int maxEpicenter=0;
				int loc = 0;
				for(int i = 0; i < ofc.Np; i++){
					if (ofc.epicenterSize[i] > maxEpicenter){
						maxEpicenter = ofc.epicenterSize[i];
						loc = i;
					}
				}
				cgMaxSizeActBoxTempAcc.accum(ofc.cg_time, maxEpicenter);
				cgMaxSizeActLocationTempAcc.accum(ofc.cg_time, loc);
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
				
				//Accums for breakdown of met calc
//				amSqAveAccum.accum(ofc.cg_time, ofc.am_sq_ave);
//				amAveSqAccum.accum(ofc.cg_time, ofc.am_ave_sq);
//				samSqAveAccum.accum(ofc.cg_time, ofc.sam_sq_ave);
//				samAveSqAccum.accum(ofc.cg_time, ofc.sam_ave_sq);

			}
			
//			System.out.println("no data pts = " + dataPointsCount + " cg time = " + ofc.cg_time + " next accum time = " + nextAccumTime);
			
			if(dataPointsCount > dataPointsPerWrite -1){ //write accumulated data
				FileUtil.initFile(sizeHistFile, params, "avalanch size histogram");
				FileUtil.printHistToFile(sizeHistFile, sizeHist);
				printAccumsToFile(iMetFile, invStressMetTempAcc, cgInvStressMetTempAcc, cgInvActivityMetTempAcc, invSizeActMetTempAcc);
				
				printAccumsToFile(sizeFile, maxSizeTempAcc,cgMaxSizeActBoxTempAcc, timeAveForMaxCGSizeActTempAcc, cgMaxSizeActLocationTempAcc);
			
				//restart counts
				maxSize = 0;
				dataPointsCount = 0;
				
				//clear temp accumulators
				invStressMetTempAcc.clear();
				cgInvStressMetTempAcc.clear();
				cgInvActivityMetTempAcc.clear();
				maxSizeTempAcc.clear();
				invSizeActMetTempAcc.clear();
				cgMaxSizeActBoxTempAcc.clear();
				timeAveForMaxCGSizeActTempAcc.clear();

			}

		
			if(ofc.plateUpdates > maxTime) Job.signalStop();
//			if(ofc.plateUpdates == 500000) Job.signalStop();
//			if(ofc.plateUpdates == 500001) Job.signalStop();
				

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
		FileUtil.initFile(sizeFile, params, "time (plate updates, fixed),  max avalanche size, max size act, time av of max size act, loc of max size act");
		sizeHistFile = params.sget("Data Dir") + File.separator + "sh.txt";	//to record size histogram data

		
		
	}

}

package rachele.damage2D.apps;

import java.awt.Color;
import java.io.File;

import rachele.damage2D.OFC_Lattice;
import rachele.util.FileUtil;
import scikit.dataset.Accumulator;
import scikit.dataset.Histogram;
import scikit.graphics.dim2.Geom2D;
import scikit.graphics.dim2.Grid;
import scikit.graphics.dim2.Plot;
import scikit.jobs.Control;
import scikit.jobs.Job;
import scikit.jobs.Simulation;
import scikit.jobs.params.DirectoryValue;

/**
* Same as basic app, 
* but without writing anything to file and with accumulators to 
* watch metrics in real time.
* 
*/
public class OFC_NoAutoWriteApp extends Simulation{

	int cg_dt;
	
	Grid grid = new Grid("Stress");
	Grid stressTimeAveGrid = new Grid ("Stress Time Ave");
	Grid cgStressGrid = new Grid("CG stress grid");
	Grid cgStressTimeAveGrid = new Grid("CG time ave stress grid");
	Grid cgGrid = new Grid("CG grid");
	Grid cgGridTimeAverage = new Grid("Time ave CG grid");
	Grid plateUpdateGrid = new Grid("Plate Update grid");
	Grid cgFailCountGrid = new Grid("cg fail count grid");
	Grid cgSizeActTimeAveGrid = new Grid("CG size act time ave");
	Grid oneEpicenter = new Grid("Special Epicienter avalanches");
	
	Plot stressMetPlot = new Plot("Stress Met Plot");
	Plot actPlot = new Plot("CG Activity I Metric");
	Plot sizePlot = new Plot("Size Plot");
	Plot cgStressMetPlot = new Plot("CG Stress Met Plot");
	Plot cgSizeActPlot = new Plot("CG Size Act Met");
	Plot sizeSumPlot = new Plot("Size sum plot");
	Plot sizeActMaxPlot = new Plot("Size act max");
	
	Plot histAll = new Plot("Histogram all sizes");
	Plot histSp = new Plot("Histogram special");
	
	OFC_Lattice ofc;
	Accumulator cgMetricAcc = new Accumulator();
	Accumulator CGstressMetAcc = new Accumulator();
	Accumulator maxSizeAcc = new Accumulator();
	Accumulator stressMetAcc = new Accumulator();
	Accumulator cgActSizeMetAcc = new Accumulator();
	Accumulator cgMaxSizeActBoxAcc = new Accumulator();
	Accumulator sizeSumAcc = new Accumulator();
	Accumulator timeAveForMaxCGSizeActAcc = new Accumulator();
	Histogram sizeHist = new Histogram(1);
	Histogram spSizeHist = new Histogram(1);
	int maxSize, fileNo;
	
	int [] avalancheOneEpicenter;
	
	public static void main(String[] args) {
		new Control(new OFC_NoAutoWriteApp(), "OFC Model");
	}
	
	public void load(Control c) {
//		c.frameTogether("Grids", stressMetPlot, cgSizeActPlot, cgGrid, cgStressMetPlot, sizePlot, cgSizeActPlot, actPlot, sizeSumPlot, sizeActMaxPlot);
		c.frameTogether("avs",plateUpdateGrid, oneEpicenter);
		c.frameTogether("size hists", histAll, histSp);
		params.add("Data Dir",new DirectoryValue("/Users/erdomi/data/damage/testRuns"));
		params.addm("Random Seed", 1);
		params.addm("CG size", 8);
		params.addm("dx", 8);
		params.addm("Coarse Grained dt (PU)", 25);
		params.addm("Equilibration Updates", 100);
		params.addm("Max Time", 1000000);
		params.addm("R", 4);// 0 -> fully connected
		params.addm("Residual Stress", 0.625);
		params.addm("Dissipation Param", 0.1);
		params.addm("Res. Max Noise", 0.125);
		params.addm("Lower Cutoff", 1);
		params.add("L");
		params.add("Time");
		params.add("Plate Updates");
		flags.add("Clear");
		flags.add("Write");
		
	}
	
	public void animate() {
		grid.setScale(0.0,1.0);
		grid.registerData(ofc.L, ofc.L, ofc.stress);
		stressTimeAveGrid.setScale(0.0,1.0);
		stressTimeAveGrid.registerData(ofc.L, ofc.L, ofc.stressTimeAve);
		cgGrid.registerData(ofc.Lp, ofc.Lp, ofc.epicenterCount);
		cgStressGrid.setScale(0.0,1.0);
		cgStressGrid.registerData(ofc.Lp, ofc.Lp, ofc.CG_Stress);
		cgGridTimeAverage.registerData(ofc.Lp, ofc.Lp, ofc.CG_ActivityTimeAve);
		plateUpdateGrid.setScale(0.0, 3);
		plateUpdateGrid.registerData(ofc.L, ofc.L, ofc.plateUpdateFailLocations);
		cgStressTimeAveGrid.setScale(0.0,1.0);
		cgStressTimeAveGrid.registerData(ofc.Lp, ofc.Lp, ofc.CG_StressTimeAve);
		cgFailCountGrid.registerData(ofc.Lp, ofc.Lp, ofc.epicenterSize);
		cgSizeActTimeAveGrid.registerData(ofc.Lp, ofc.Lp, ofc.CG_SizeActTimeAve);
		
		grid.clearDrawables();
		double radius = 1.0/(2.0*ofc.L);
		double failSite_y = ((double)(ofc.epicenterSite/ofc.L))/ofc.L + radius;
		double failSite_x = ((double)(ofc.epicenterSite%ofc.L))/ofc.L + radius;
		grid.addDrawable(
				Geom2D.circle(failSite_x, failSite_y, radius, Color.GREEN));
		
		oneEpicenter.registerData(ofc.L, ofc.L, avalancheOneEpicenter);
		
		plateUpdateGrid.clearDrawables();
		plateUpdateGrid.addDrawable(
				Geom2D.circle(failSite_x, failSite_y, radius, Color.GREEN));
		
		actPlot.setAutoScale(true);
		actPlot.registerLines("Activity Metric", cgMetricAcc, Color.BLACK);
		sizePlot.setAutoScale(true);
		sizePlot.registerPoints("Maz EQ size", maxSizeAcc, Color.RED);
		sizePlot.registerPoints("av size sum", sizeSumAcc, Color.BLUE);
		sizePlot.registerPoints("Max Size Act Block", cgMaxSizeActBoxAcc, Color.pink);
		cgSizeActPlot.setAutoScale(true);
		cgSizeActPlot.registerLines("size act", cgActSizeMetAcc, Color.green);
		cgStressMetPlot.setAutoScale(true);
		cgStressMetPlot.registerLines("CG stress met", CGstressMetAcc, Color.BLUE);
		stressMetPlot.setAutoScale(true);
		stressMetPlot.registerLines("stress met", stressMetAcc, Color.BLACK);
		sizeSumPlot.setAutoScale(true);
		sizeSumPlot.registerPoints("av size sum", sizeSumAcc, Color.BLUE);
		sizeActMaxPlot.setAutoScale(true);
		sizeActMaxPlot.registerPoints("Max Size Act Block", cgMaxSizeActBoxAcc, Color.pink);
		sizeActMaxPlot.registerPoints("Max Size Act Block Time Ave", timeAveForMaxCGSizeActAcc, Color.black);
		
		histAll.registerBars("all sizes", sizeHist, Color.BLACK);
		histSp.registerBars("Special sizes", spSizeHist, Color.RED);
		
		params.set("Time", ofc.cg_time);
		params.set("Plate Updates", ofc.plateUpdates);
		
		if(flags.contains("Write")){
			writeCgStressTimeAveGrid();
		}
		if (flags.contains("Clear")){
			cgMetricAcc.clear();
			maxSizeAcc.clear();
			CGstressMetAcc.clear();
			stressMetAcc.clear();
			cgActSizeMetAcc.clear();
			cgMaxSizeActBoxAcc.clear();
			sizeSumAcc.clear();
			timeAveForMaxCGSizeActAcc.clear();
		}
//		flags.clear();
	}

	public void clear() {
	}

	public void run() {
		ofc = new OFC_Lattice(params);
		cg_dt = params.iget("Coarse Grained dt (PU)");
		double nextRecordTime = 0;
		int maxTime = params.iget("Max Time");
		int sizeSum=0;
		maxSize = 0;
		fileNo = 0;
		int spAvNo = 0;
		
		avalancheOneEpicenter = new int [ofc.N];
		int specialEpicenter = ofc.N/2 + ofc.L/2;
		
		//equilibrate
		ofc.initEquilibrate(params.iget("Equilibration Updates"));

		while (ofc.plateUpdates < 0){
			ofc.equilibrate();
			Job.animate();
			if(ofc.plateUpdates == 0) Job.signalStop();
		}
		
		while(true){

			ofc.step();
			ofc.cgFailCountAdd();
			

			
			int size = ofc.avSize;
			sizeHist.accum(size);
			if (size > maxSize) {
				maxSize = size;
			}
			sizeSum += size;
			if (ofc.epicenterSite == specialEpicenter){
				spAvNo += 1;
				System.out.println("Av " + spAvNo + " at time " + ofc.plateUpdates );
				spSizeHist.accum(size);
				for (int i = 0; i < ofc.N; i++){
					avalancheOneEpicenter[i] += ofc.plateUpdateFailLocations[i];
				}
			}
			
			double iMet = 1.0/ofc.calcStressMetric();

			if(ofc.plateUpdates > nextRecordTime){
				
				int maxEpicenter=0;
				int loc = -1;
				for(int i = 0; i < ofc.Np; i++){
					if (ofc.epicenterSize[i] > maxEpicenter){
						maxEpicenter = ofc.epicenterSize[i];
						loc = i;
					}
				}
				cgMaxSizeActBoxAcc.accum(ofc.cg_time, maxEpicenter);
				timeAveForMaxCGSizeActAcc.accum(ofc.cg_time, ofc.CG_SizeActTimeAve[loc]);
				
				//stress metric
				stressMetAcc.accum(ofc.cg_time, iMet);
				
				//CG stress metric
				double cgInverseStressMetric = 1.0/ofc.calcCG_stressMetric();
				CGstressMetAcc.accum(ofc.cg_time, cgInverseStressMetric);
				
				
				double activityOmega = ofc.calcCG_activityMetric();

				double cgInverseActivityMetric = 1.0/activityOmega;
				cgMetricAcc.accum(ofc.cg_time, cgInverseActivityMetric);

				double cgInvSizeActMet = 1.0/ofc.calcCG_sizeActMetric();
				cgActSizeMetAcc.accum(ofc.cg_time, cgInvSizeActMet);
				
				maxSizeAcc.accum(ofc.cg_time, maxSize);

				sizeSumAcc.accum(ofc.cg_time, sizeSum);
				
				nextRecordTime += cg_dt;
				maxSize = 0;
				sizeSum = 0;
			}


			if(ofc.cg_time > maxTime) Job.signalStop();
				
			Job.animate();
		}
	}
	
	void writeCgStressTimeAveGrid(){
		String fileName = params.sget("Data Dir") + File.separator + "g" + fileNo;  
		FileUtil.initFile(fileName, params);
		for (int i=0; i < ofc.Np; i++){
			int x = i%ofc.Lp;
			int y = i/ofc.Lp;
			FileUtil.printlnToFile(fileName, x, y,ofc.CG_ActivityTimeAve[i]);
		}
		fileNo +=1;
	}
	
}

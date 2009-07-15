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
	Grid cgGrid = new Grid(" CG grid");
	Grid cgGridTimeAverage = new Grid("Time ave CG grid");
	Grid plateUpdateGrid = new Grid("Plate Update grid");
	Plot stressMetPlot = new Plot("Stress Met Plot");
	Plot actPlot = new Plot("CG Activity I Metric");
	Plot sizePlot = new Plot("Size Plot");
	Plot cgStressMetPlot = new Plot("CG Stress Met Plot");
	OFC_Lattice ofc;
	Accumulator cgMetricAcc = new Accumulator();
	Accumulator CGstressMetAcc = new Accumulator();
	Accumulator maxSizeAcc = new Accumulator();
	Accumulator stressMetAcc = new Accumulator();
	Histogram sizeHist = new Histogram(1);
	int maxSize, fileNo;
	
	public static void main(String[] args) {
		new Control(new OFC_NoAutoWriteApp(), "OFC Model");
	}
	
	public void load(Control c) {
		c.frameTogether("Grids", grid, stressTimeAveGrid, plateUpdateGrid, cgStressGrid, cgStressTimeAveGrid, sizePlot, cgStressMetPlot, actPlot, stressMetPlot);
		params.add("Data Dir",new DirectoryValue("/home/erdomi/data/damage/testRuns"));
		params.addm("Random Seed", 1);
		params.addm("CG size", 8);
		params.addm("dx", 8);
		params.addm("Coarse Grained dt", 1);
		params.addm("Equilibration Updates", 1000);
		params.addm("Max Time", 1000000);
		params.addm("R", 0);// 0 -> fully connected
		params.addm("Residual Stress", 0.625);
		params.addm("Dissipation Param", 0.05);
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
		
		grid.clearDrawables();
		double radius = 1.0/(2.0*ofc.L);
		double failSite_y = ((double)(ofc.epicenterSite/ofc.L))/ofc.L + radius;
		double failSite_x = ((double)(ofc.epicenterSite%ofc.L))/ofc.L + radius;
		grid.addDrawable(
				Geom2D.circle(failSite_x, failSite_y, radius, Color.GREEN));
		
		plateUpdateGrid.clearDrawables();
		plateUpdateGrid.addDrawable(
				Geom2D.circle(failSite_x, failSite_y, radius, Color.GREEN));
		
		actPlot.setAutoScale(true);
		actPlot.registerLines("Activity Metric", cgMetricAcc, Color.BLACK);
		sizePlot.setAutoScale(true);
		sizePlot.registerPoints("Maz EQ size", maxSizeAcc, Color.RED);
		cgStressMetPlot.setAutoScale(true);
		cgStressMetPlot.registerLines("CG stress met", CGstressMetAcc, Color.BLUE);
		stressMetPlot.setAutoScale(true);
		stressMetPlot.registerLines("stress met", stressMetAcc, Color.BLACK);
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
		}
		flags.clear();
	}

	public void clear() {
	}

	public void run() {
		ofc = new OFC_Lattice(params);
		cg_dt = params.iget("Coarse Grained dt");
		double nextRecordTime = 0;
		int maxTime = params.iget("Max Time");
		maxSize = 0;
		fileNo = 0;
		
		//equilibrate
		ofc.initEquilibrate(params.iget("Equilibration Updates"));

		while (ofc.plateUpdates < 0){
			ofc.equilibrate();
//			System.out.println("init done, R = " + ofc.R);
			Job.animate();
			if(ofc.plateUpdates == 0) Job.signalStop();
		}
		
		while(true){

			ofc.step();

			int size = ofc.avSize;
			sizeHist.accum(size);
			if (size > maxSize) {
				maxSize = size;
			}
//			double iMet = ofc.calcInverseMetric();
			double iMet = 1.0/ofc.calcStressMetric();

			if(ofc.plateUpdates > nextRecordTime){
//				System.out.println(ofc.plateUpdates + " PU " + "nextRecordTime = " + nextRecordTime);
				
				//stress metric
				stressMetAcc.accum(ofc.cg_time, iMet);
				
				//CG stress metric
				double cgInverseStressMetric = 1.0/ofc.calcCG_stressMetric();
				CGstressMetAcc.accum(ofc.cg_time, cgInverseStressMetric);
				
				
				double activityOmega = ofc.calcCG_activityMetric();

				double cgInverseActivityMetric = 1.0/activityOmega;
				cgMetricAcc.accum(ofc.cg_time, cgInverseActivityMetric);

				maxSizeAcc.accum(ofc.cg_time, maxSize);

				nextRecordTime += cg_dt;
				maxSize = 0;
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

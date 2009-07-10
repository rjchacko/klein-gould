package rachele.damage2D.apps;

import java.awt.Color;
import rachele.damage2D.OFC_Lattice;
import scikit.dataset.Accumulator;
import scikit.dataset.Histogram;
import scikit.graphics.dim2.Geom2D;
import scikit.graphics.dim2.Grid;
import scikit.graphics.dim2.Plot;
import scikit.jobs.Control;
import scikit.jobs.Job;
import scikit.jobs.Simulation;

/**
* Same as basic app, 
* but without writing anything to file and with accumulators to 
* watch metrics in real time.
* 
*/
public class OFC_AppNoFileWrite extends Simulation{

	int cg_dt;
	
	Grid grid = new Grid("Lattice");
	Grid cgGrid = new Grid(" CG grid");
	Grid cgGridTimeAverage = new Grid("Time ave CG grid");
	Grid plateUpdateGrid = new Grid("Plate Update grid");
//	Plot sizePlot = new Plot("Size Histogram");
	Plot actPlot = new Plot("CG Activity I Metric");
	Plot sizePlot = new Plot("Size Plot");
	Plot stressMetPlot = new Plot("Stress Met Plot");
	OFC_Lattice ofc;
	Accumulator cgMetricAcc = new Accumulator();
	Accumulator maxSizeAcc = new Accumulator();
	Accumulator CGstressMetAcc = new Accumulator();
	Accumulator sizeAcc = new Accumulator();
	Histogram sizeHist = new Histogram(1);
	String iMetFile;
//	String sizeFile;
	String sizeHistFile;
	int maxSize;
//	String cgInvMetricFile;
	
	public static void main(String[] args) {
		new Control(new OFC_AppNoFileWrite(), "OFC Model");
	}
	
	public void load(Control c) {
		c.frameTogether("Grids", grid, cgGrid, cgGridTimeAverage, plateUpdateGrid, sizePlot, stressMetPlot, actPlot);
//		c.frameTogether("Data", sizePlot, stressMetPlot, actPlot);
//		params.add("Data Dir",new DirectoryValue("/Users/erdomi/data/damage/testRuns"));
		params.addm("Random Seed", 1);
		params.addm("CG size", 32);
		params.addm("dx", 1);
		params.addm("Coarse Grained dt", 50);
		params.addm("Equilibration Updates", 1000);
		params.addm("Max Time", 1000000);
		params.addm("R", 16);// 0 -> fully connected
		params.addm("Residual Stress", 0.625);
		params.addm("Dissipation Param", 0.05);
		params.addm("Res. Max Noise", 0.125);
		params.addm("Lower Cutoff", 1);
		params.add("L");
		params.add("Time");
		params.add("Plate Updates");
		flags.add("Clear");
		
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
		
		plateUpdateGrid.clearDrawables();
		plateUpdateGrid.addDrawable(
				Geom2D.circle(failSite_x, failSite_y, radius, Color.GREEN));
		
		actPlot.setAutoScale(true);
		actPlot.registerLines("Activity Metric", cgMetricAcc, Color.BLACK);
		sizePlot.setAutoScale(true);
		sizePlot.registerPoints("Maz EQ size", sizeAcc, Color.RED);
		stressMetPlot.setAutoScale(true);
		stressMetPlot.registerLines("stress met", CGstressMetAcc, Color.BLUE);
//		iterPlot.setAutoScale(true);
//		iterPlot.registerLines("Iterations", sizeAcc, Color.RED);
		
//		metricPlot.setAutoScale(true);
//		metricPlot.registerLines("Metric", inverseMetricAcc, Color.BLACK);
//		sizePlot.setAutoScale(true);
//		sizePlot.registerBars("Size", sizeHist, Color.BLACK);
		params.set("Time", ofc.time);
		params.set("Plate Updates", ofc.plateUpdates);
		
		if (flags.contains("Clear")){
			cgMetricAcc.clear();
			sizeAcc.clear();
			CGstressMetAcc.clear();
			flags.clear();
		}
	}

	public void clear() {
	}

	public void run() {
		ofc = new OFC_Lattice(params);
//		initFiles();
		cg_dt = params.iget("Coarse Grained dt");
//		boolean prestep = true;
		double nextRecordTime = 0;
//		boolean findActivityOmega0 = true; 
		int maxTime = params.iget("Max Time");
		maxSize = 0;
//		double activityOmega0=0;
//		double CGstressOmega0 = 0;
		
		//equilibrate
		ofc.initEquilibrate(params.iget("Equilibration Updates"));
		while (ofc.plateUpdates < 0){
			ofc.equilibrate();
			Job.animate();
			if(ofc.plateUpdates == 0) Job.signalStop();
		}
		
		while(true){
//			if(prestep){
//			ofc.prestep();
//			prestep =false;
//			}else{	
			ofc.clearPlateUpdateFileLocs();
			ofc.step();
//			prestep = true;

			int size = ofc.avSize;
			sizeHist.accum(size);
//			sizeStore.accum(ofc.time, size);
			if (size > maxSize) {
				maxSize = size;
//				System.out.print("max size " + maxSize);
			}
			double iMet = ofc.calcInverseMetric();

			if(ofc.time > nextRecordTime){

//				FileUtil.printlnToFile(sizeFile, ofc.time, size);
//				FileUtil.initFile(sizeHistFile, params, "avalanch size histogram");
//				FileUtil.printHistToFile(sizeHistFile, sizeHist);
				double activityOmega = ofc.calcCG_activityMetric();
				double CGstressOmega = ofc.calcCG_stressMetric();
//				if(findActivityOmega0){
//					if(activityOmega != 0){
//						activityOmega0 = activityOmega;
//						findActivityOmega0 = false;
//					}
//					System.out.println("omega0 = " + activityOmega0 + " at time = " + ofc.time);
//					FileUtil.printlnToFile(iMetFile, "# Omega(t=0) = ", activityOmega0);
//				}
//				if(nextRecordTime == 0){
//					CGstressOmega0 = CGstressOmega;
//					System.out.println("omega0  stress= " + CGstressOmega0);
//				}

				double reducedTime = ofc.time/cg_dt;
				double cgInverseActivityMetric = 1.0/activityOmega;
				cgMetricAcc.accum(reducedTime, cgInverseActivityMetric);
				double cgInverseStressMetric = 1.0/CGstressOmega;
				CGstressMetAcc.accum(reducedTime, cgInverseStressMetric);
				sizeAcc.accum(reducedTime, maxSize);
//				FileUtil.printlnToFile(iMetFile, ofc.time, iMet, reducedTime, cgInverseActivityMetric, cgInverseStressMetric, maxSize);
//				FileUtil.printAccumToFileNoErrorBars(sizeFile, sizeStore);
//				FileUtil.printlnToFile(cgInvMetricFile, ofc.time, cgInverseMetric);
				nextRecordTime += cg_dt;
				maxSize = 0;
//				sizeStore.clear();
			}
//			}


			if(ofc.time > maxTime) Job.signalStop();
				
			Job.animate();
		}
	}
	
//	void initFiles(){
//		iMetFile = params.sget("Data Dir") + File.separator + "im.txt";  // to record iverse metric data
//		FileUtil.initFile(iMetFile, params, " time (plate updates), stress inverse metric, time/coarse grained time, coarse grained activity metric, coarse grained stress metric, size of avalanche");
////		sizeFile = params.sget("Data Dir") + File.separator + "s.txt";   //to record size vs time data
////		FileUtil.initFile(sizeFile, params, "avalanch size vs time");
//		sizeHistFile = params.sget("Data Dir") + File.separator + "sh.txt";	//to record size histogram data
////		cgInvMetricFile = params.sget("Data Dir") + File.separator + "cg.txt";
////		FileUtil.initFile(cgInvMetricFile, params, "Coarse Grained Inverse Metric");
//	}

}

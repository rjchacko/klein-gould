package rachele.damage2D.apps;

import java.awt.Color;
import java.io.File;

import rachele.damage2D.OFC_Lattice;
import rachele.util.FileUtil;
//import scikit.dataset.Accumulator;
import scikit.dataset.Histogram;
import scikit.graphics.dim2.Geom2D;
import scikit.graphics.dim2.Grid;
//import scikit.graphics.dim2.Plot;
import scikit.jobs.Control;
import scikit.jobs.Job;
import scikit.jobs.Simulation;
import scikit.jobs.params.DirectoryValue;
import scikit.util.Utilities;

/**
* Basic app.  
* 
* Includes event size lower cutoff option.  
* 
* Outputs every CG time chunk (=cg_dt) stress metric, coarse grained stress metric, 
* coarse grained activity metric, max size of cg_dt time chunk and a (GR) size histogram
* 
* To run damage mode, use OFC_DamageApp.
*/
public class OFC_App extends Simulation{

	int cg_dt;
	
	Grid grid = new Grid("Lattice");
	Grid cgGrid = new Grid(" CG grid");
	Grid cgGridTimeAverage = new Grid("Time ave CG grid");
//	Plot sizePlot = new Plot("Size Histogram");
	OFC_Lattice ofc;
//	Accumulator cgMetricAcc = new Accumulator();
//	Accumulator sizeStore = new Accumulator();  //To store data at each plate update within a bracket
	
	Histogram sizeHist = new Histogram(1);
	String iMetFile;
//	String sizeFile;
	String sizeHistFile;
	int maxSize;
//	String cgInvMetricFile;
	
	public static void main(String[] args) {
		new Control(new OFC_App(), "OFC Model");
	}
	
	public void load(Control c) {
		c.frameTogether("Grids", grid, cgGrid, cgGridTimeAverage);
//		c.frameTogether("Data", sizePlot);
		params.add("Data Dir",new DirectoryValue("/Users/erdomi/data/damage/testRuns"));
		params.addm("Random Seed", 1);
		params.addm("CG size", 32);
		params.addm("dx", 8);
		params.addm("Coarse Grained dt", 500);
		params.addm("Equilibration Updates", 1000000);
		params.addm("Max Time", 1000000);
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
		cgGridTimeAverage.registerData(ofc.Lp, ofc.Lp, ofc.CG_ActivityTimeAve);
		
		grid.clearDrawables();
		double radius = 1.0/(2.0*ofc.L);
		double failSite_y = ((double)(ofc.epicenterSite/ofc.L))/ofc.L + radius;
		double failSite_x = ((double)(ofc.epicenterSite%ofc.L))/ofc.L + radius;
		grid.addDrawable(
				Geom2D.circle(failSite_x, failSite_y, radius, Color.GREEN));
		
//		iterPlot.setAutoScale(true);
//		iterPlot.registerLines("Iterations", sizeAcc, Color.RED);
		
//		metricPlot.setAutoScale(true);
//		metricPlot.registerLines("Metric", inverseMetricAcc, Color.BLACK);
//		sizePlot.setAutoScale(true);
//		sizePlot.registerBars("Size", sizeHist, Color.BLACK);
		params.set("Time", Utilities.format(ofc.cg_time));
		params.set("Plate Updates", ofc.plateUpdates);
	}

	public void clear() {
	}

	public void run() {
		ofc = new OFC_Lattice(params);
		initFiles();
		cg_dt = params.iget("Coarse Grained dt");
//		boolean prestep = true;
		double nextRecordTime = 0;
		boolean findActivityOmega0 = true; 
		int maxTime = params.iget("Max Time");
		maxSize = 0;
		double activityOmega0=0;
		double CGstressOmega0 = 0;
		
		//equilibrate
		ofc.initEquilibrate(params.iget("Equilibration Updates"));
		while (ofc.plateUpdates < 0){
			ofc.equilibrate();
			Job.animate();
		}
		
		while(true){
//			if(prestep){
//			ofc.prestep();
//			prestep =false;
//			}else{				
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

			if(ofc.cg_time > nextRecordTime){

//				FileUtil.printlnToFile(sizeFile, ofc.time, size);
				FileUtil.initFile(sizeHistFile, params, "avalanch size histogram");
				FileUtil.printHistToFile(sizeHistFile, sizeHist);
				double activityOmega = ofc.calcCG_activityMetric();
				double CGstressOmega = ofc.calcCG_stressMetric();
				if(findActivityOmega0){
					if(activityOmega != 0){
						activityOmega0 = activityOmega;
						findActivityOmega0 = false;
					}
					System.out.println("omega0 = " + activityOmega0 + " at time = " + ofc.cg_time);
					FileUtil.printlnToFile(iMetFile, "# Omega(t=0) = ", activityOmega0);
				}
				if(nextRecordTime == 0){
					CGstressOmega0 = CGstressOmega;
					System.out.println("omega0  stress= " + CGstressOmega0);
				}
				double cgInverseActivityMetric = activityOmega0/activityOmega;
				double cgInverseStressMetric = CGstressOmega0/CGstressOmega;
				double reducedTime = ofc.cg_time/cg_dt;
				FileUtil.printlnToFile(iMetFile, ofc.cg_time, iMet, reducedTime, cgInverseActivityMetric, cgInverseStressMetric, maxSize);
//				FileUtil.printAccumToFileNoErrorBars(sizeFile, sizeStore);
//				FileUtil.printlnToFile(cgInvMetricFile, ofc.time, cgInverseMetric);
				nextRecordTime += cg_dt;
				maxSize = 0;
//				sizeStore.clear();
			}
//			}

			if(ofc.cg_time > maxTime) Job.signalStop();
				
			Job.animate();
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

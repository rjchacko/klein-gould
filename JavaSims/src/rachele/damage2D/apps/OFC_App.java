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
import scikit.util.Utilities;

public class OFC_App extends Simulation{

	double cg_dt;
	
	Grid grid = new Grid("Lattice");
	Grid cgGrid = new Grid(" CG grid");
//	Plot iterPlot = new Plot("Iterations");
//	Plot metricPlot = new Plot("Metric Plot");
	Plot sizePlot = new Plot("Size Histogram");
	OFC_Lattice ofc;
	Accumulator cgMetricAcc;
//	Accumulator sizeAcc = new Accumulator(100);
	Histogram sizeHist = new Histogram(1);
//	Accumulator inverseMetricAcc = new Accumulator();
	String iMetFile;
	String sizeFile;
	String sizeHistFile;
	
	public static void main(String[] args) {
		new Control(new OFC_App(), "OFC Model");
	}
	
	public void load(Control c) {
		c.frameTogether("Grids", grid, cgGrid);
		c.frameTogether("Data", sizePlot);
		params.add("Data Dir",new DirectoryValue("/Users/erdomi/data/damage/testRuns"));
		params.addm("Random Seed", 1);
		params.addm("CG size", 256);
		params.addm("dx", 1);
		params.addm("Coarse Grained dt", 100);
		params.addm("Equilibration Updates", 10000);
		params.addm("R", 0);// 0 -> fully connected
		params.addm("Residual Stress", 0.625);
		params.addm("Dissipation Param", 0.2);
		params.addm("Res. Max Noise", 0.125);
		params.add("L");
		params.add("Time");
		params.add("Plate Updates");
	}
	
	public void animate() {
		grid.registerData(ofc.L, ofc.L, ofc.stress);
		cgGrid.registerData(ofc.Lp, ofc.Lp, ofc.cgCount);
		
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
		sizePlot.setAutoScale(true);
		sizePlot.registerBars("Size", sizeHist, Color.BLACK);
		params.set("Time", Utilities.format(ofc.time));
		params.set("Plate Updates", ofc.plateUpdates);
	}

	public void clear() {
	}

	public void run() {
		ofc = new OFC_Lattice(params);
		initFiles();
		cg_dt = params.fget("Coarse Grained dt");
		cgMetricAcc = new Accumulator();
		boolean prestep = true;
		double nextRecordTime = cg_dt;
		
		//equilibrate
		ofc.initEquilibrate(params.iget("Equilibration Updates"));
		while (ofc.plateUpdates < 0){
			ofc.equilibrate();
			Job.animate();
		}
		
		while(true){
			if(prestep){
				ofc.prestep();
				prestep =false;
			}else{				
				ofc.step();
				prestep = true;

				int size = ofc.avSize;

				FileUtil.printlnToFile(iMetFile, ofc.time, ofc.calcInverseMetric());
				FileUtil.printlnToFile(sizeFile, ofc.time, size);
				sizeHist.accum(size);

				if(ofc.time > nextRecordTime){
					FileUtil.initFile(sizeHistFile, params, "avalanch size histogram");
					FileUtil.printHistToFile(sizeHistFile, sizeHist);
//					double cgMetric = ofc.calcCG_Metric();
					nextRecordTime += cg_dt;
				}
			}


			Job.animate();
		}
	}
	
	void initFiles(){
		iMetFile = params.sget("Data Dir") + File.separator + "im.txt";  // to record iverse metric data
		String message = "inverse metric vs time";
		FileUtil.initFile(iMetFile, params, message);
		sizeFile = params.sget("Data Dir") + File.separator + "s.txt";   //to record size vs time data
		message = "avalanch size vs time";
		FileUtil.initFile(sizeFile, params, message);
		sizeHistFile = params.sget("Data Dir") + File.separator + "sh.txt";	//to record size histogram data
	}

}

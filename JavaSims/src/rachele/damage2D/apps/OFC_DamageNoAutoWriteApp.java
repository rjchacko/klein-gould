package rachele.damage2D.apps;

import java.awt.Color;
//import java.text.Format;
//import java.io.File;
import rachele.damage2D.OFC_DamageLattice;
//import rachele.util.FileUtil;
import scikit.dataset.Accumulator;
import scikit.dataset.Histogram;
import scikit.graphics.dim2.Geom2D;
import scikit.graphics.dim2.Grid;
import scikit.jobs.Control;
import scikit.jobs.Job;
import scikit.jobs.Simulation;
import scikit.jobs.params.ChoiceValue;
//import scikit.jobs.params.DirectoryValue;
import scikit.util.Utilities;

public class OFC_DamageNoAutoWriteApp extends Simulation{

	int cg_dt;

	Grid grid = new Grid("Lattice");
	Grid cgGrid = new Grid("CG grid");
	Grid deadGrid = new Grid("Alive/Dead Lattice");
	//	Plot sizePlot = new Plot("Size Histogram");
	OFC_DamageLattice ofc;
	Accumulator cgMetricAcc = new Accumulator();
	//	Accumulator sizeStore = new Accumulator();  //To store data at each plate update within a bracket

	Histogram sizeHist = new Histogram(1);
	String iMetFile;
	//	String sizeFile;
	String sizeHistFile;
	int maxSize;
	//	String cgInvMetricFile;

	public static void main(String[] args) {
		new Control(new OFC_DamageNoAutoWriteApp(), "OFC Damage Model No Write");
	}
	

	public void load(Control c) {
		c.frameTogether("Grids", grid, cgGrid, deadGrid);
		params.add("Interaction", new ChoiceValue("Circle", "Fully Connected", "Square", "Small World") );
		params.addm("Random Seed", 1);
		params.addm("CG size", 30);
		params.addm("dx", 9);
		params.addm("Coarse Grained dt", 100);
		params.addm("Equilibration Updates", 1000);
		params.addm("R", 10);// 0 -> fully connected
		params.addm("Residual Stress", 0.625);
		params.addm("Dissipation Param", 0.3);
		params.addm("Res. Max Noise", 0.125);
		params.addm("Lower Cutoff", 1);
		params.addm("Mean Max Failures", 1);
		params.addm("Failures Max Noise", 0);
		params.addm("Mean Heal Time", 1);
		params.addm("Heal Time Noise", 0);
		params.add("L");
		params.add("Time");
		params.add("Av Size");
		params.add("Plate Updates");
		params.add("Percent dead sites");
	}

	public void animate() {
		grid.setScale(0.0, 1.0);
		grid.registerData(ofc.L, ofc.L, ofc.stress);
		cgGrid.registerData(ofc.Lp, ofc.Lp, ofc.epicenterCount);
		int maxNoFails = params.iget("Mean Max Failures") + params.iget("Failures Max Noise");
		deadGrid.setScale(-1, maxNoFails-1);
		deadGrid.registerData(ofc.L, ofc.L, ofc.noFails);
		
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
		double percentDead = ofc.noDeadSites/(double)(ofc.L*ofc.L);
		params.set("Percent dead sites",  Utilities.format(percentDead));
		params.set("Av Size", ofc.avSize);
		if (percentDead >= 1.0) Job.signalStop();
	}

	public void clear() {
	}

	
	public void run() {
		ofc = new OFC_DamageLattice(params, "damage");
		cg_dt = params.iget("Coarse Grained dt");
		maxSize = 0;

		//equilibrate
		ofc.initEquilibrate(params.iget("Equilibration Updates"));
		while (ofc.plateUpdates < 0){
			ofc.equilibrate();
			Job.animate();
		}

		while(true){
			ofc.healPreStep();
			Job.animate();
			while (ofc.nextSiteToFail >= 0){
				ofc.healStepIter();
				Job.animate();
			}
			if (ofc.avSize >= ofc.lowerCutoff) 	ofc.epicenterCount[ofc.findCG_site(ofc.epicenterSite, ofc.dx)] +=1;

			int size = ofc.avSize;
//			System.out.println(" av size = " + size + " at time " + ofc.plateUpdates);
			sizeHist.accum(size);
			if (size > maxSize) {
				maxSize = size;
			}

			Job.animate();
		}
	}


}

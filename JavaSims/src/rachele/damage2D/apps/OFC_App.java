package rachele.damage2D.apps;

import java.awt.Color;

import rachele.damage2D.OFC_Lattice;
import scikit.graphics.dim2.Geom2D;
import scikit.graphics.dim2.Grid;
import scikit.graphics.dim2.Plot;
import scikit.jobs.Control;
import scikit.jobs.Job;
import scikit.jobs.Simulation;

public class OFC_App extends Simulation{

	Grid grid = new Grid("Lattice");
	Grid cgGrid = new Grid(" CG grid");
	Plot iterPlot = new Plot("Iterations");
	Plot metricPlot = new Plot("Metric Plot");
	OFC_Lattice ofc;
	
	public static void main(String[] args) {
		new Control(new OFC_App(), "OFC Model");
	}
	
	public void load(Control c) {
		c.frameTogether("Data", grid, iterPlot, cgGrid, metricPlot);
		params.addm("Random Seed", 1);
		params.addm("CG size", 512);
		params.addm("dx", 1);
		params.addm("dt", 100);
		params.addm("R", 2);
		params.addm("Residual Stress", 0.625);
		params.addm("Dissipation Param", 0.2);
		params.addm("Res. Max Noise", 0.125);
		params.add("L");
		params.add("Time");
		params.add("CG Time");
	}
	
	public void animate() {
		grid.registerData(ofc.L, ofc.L, ofc.stress);
		cgGrid.registerData(ofc.Lp, ofc.Lp, ofc.cgCount);
		
		grid.clearDrawables();
		double radius = 1.0/(2.0*ofc.L);
		double failSite_y = ((double)(ofc.failSite/ofc.L))/ofc.L + radius;
		double failSite_x = ((double)(ofc.failSite%ofc.L))/ofc.L + radius;
		grid.addDrawable(
				Geom2D.circle(failSite_x, failSite_y, radius, Color.GREEN));
		iterPlot.registerLines("Iterations", ofc.iterAcc, Color.RED);
		metricPlot.registerLines("Metric", ofc.inverseMetricAcc, Color.BLACK);
		params.set("Time", ofc.time);
	}

	public void clear() {

		
	}


	public void run() {
		ofc = new OFC_Lattice(params);
		boolean prestep = true;
		while(true){
			if(prestep){
				ofc.prestep();
				prestep =false;
			}else{				
				ofc.step();
			prestep = true;
			}
			Job.animate();
//			if(ofc.time % 1000 == 0){
//				System.out.println("time = " + ofc.time + " metric sum = " + ofc.metricSum);
//			}
		}
		
	}

}

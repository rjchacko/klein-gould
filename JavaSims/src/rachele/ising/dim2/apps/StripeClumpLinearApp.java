package rachele.ising.dim2.apps;

import java.awt.Color;

import scikit.dataset.PointSet;
import scikit.graphics.dim2.Plot;
import scikit.jobs.Control;
import scikit.jobs.Job;
import scikit.jobs.Simulation;
import scikit.jobs.params.ChoiceValue;
import scikit.jobs.params.DirectoryValue;
import scikit.jobs.params.DoubleValue;
import scikit.jobs.params.FileValue;
import rachele.ising.dim2.IsingField2D;
import rachele.ising.dim2.StripeClumpFieldSim;

public class StripeClumpLinearApp extends Simulation{
	Plot fkPlot = new Plot ("f(x) plot");
	Plot eigenvector1plot = new Plot ("EV1 Plot");
	Plot eigenvector2plot = new Plot ("EV2 Plot");
	Plot eigenvector3plot = new Plot ("EV2 Plot");
	Plot eigenvalues = new Plot ("Eigenvalues");
	Plot phiPlot = new Plot("Phi Plot");
	IsingField2D ising;
    StripeClumpFieldSim sc;
    
	public static void main(String[] args) {
		new Control(new StripeClumpLinearApp(), "Stripe -> Clump");
	}
	
	public void animate() { 
		phiPlot.registerLines("phi 0", new PointSet( 0, 1, sc.phi0), Color.black);
		fkPlot.registerLines("f_k", new PointSet( 0, 1, sc.f_k), Color.red);
		eigenvalues.registerLines("eigenvalues", new PointSet( 0, 1, sc.eigenvalue), Color.blue);
		eigenvector1plot.registerLines("eigenvec 1", new PointSet( 0, 1, sc.VV[ising.Lp-1]), Color.yellow);
		eigenvector2plot.registerLines("eigenvec 2", new PointSet( 0, 1, sc.VV[ising.Lp-2]), Color.green);
		eigenvector3plot.registerLines("eigenvec 3", new PointSet( 0, 1, sc.VV[ising.Lp-3]), Color.orange);
	}
	
	public void load(Control c) {
		params.add("1D Input File", new FileValue("/home/erdomi/data/lraim/configs1dAutoName/config"));
		params.add("Data Dir", new DirectoryValue("/home/erdomi/data/lraim/stripeToClumpInvestigation/ftResults/StripeClumpLinearApp"));
		params.addm("Interaction", new ChoiceValue("Square", "Circle"));
		params.addm("Dynamics?", new ChoiceValue("Langevin No M Convervation"));
		params.add("Init Conditions", new ChoiceValue("Read 1D Soln", "Read From File","Random Gaussian"));
		params.addm("Approx", new ChoiceValue("None", "Modified Dynamics"));
		params.addm("Noise", new DoubleValue(1.0, 0.0, 1.0).withSlider());
		params.addm("Horizontal Slice", new DoubleValue(0.5, 0, 0.9999).withSlider());
		params.addm("Vertical Slice", new DoubleValue(0.5, 0, 0.9999).withSlider());
		params.addm("kR", new DoubleValue(5.135622302, 0.0, 6.0).withSlider());
		params.addm("T", 0.04);
		params.addm("H", 0.80);
		params.addm("tolerance", 0.0001);
		params.addm("J", -1.0);
		params.addm("R", 2000000.0);
		params.addm("Random seed", 0);
		params.add("L/R", 2.7826087);
		params.add("R/dx", 50.0);
		params.add("kR bin-width", 0.1);
		params.add("Magnetization", 0.0);
		params.addm("ky", 2);
		params.addm("Max Time", 100.0);
		params.addm("dt", 0.005);
		c.frameTogether("f(x) vs phi(x)", phiPlot, fkPlot, eigenvalues, eigenvector1plot, eigenvector2plot, eigenvector3plot);
	}


	public void clear() {
	}

	public void run() {
		ising = new IsingField2D(params);
		sc = new StripeClumpFieldSim(ising, params);
		Job.animate();
	}
}
package rachele.phi4.apps;

import java.awt.Color;

import rachele.phi4.phi4;
import rachele.util.FourierTransformer;
import scikit.dataset.Accumulator;
import scikit.graphics.dim2.Grid;
import scikit.graphics.dim2.Plot;
import scikit.jobs.Control;
import scikit.jobs.Job;
import scikit.jobs.Simulation;
import scikit.jobs.params.ChoiceValue;
import scikit.jobs.params.DoubleValue;
/**
* 
* The purpose of this app is to investigate the results of a paper by 
* Corberi, Coniglio, and Zanetti
*
*/
public class Phi4App extends Simulation{

	phi4 field;
	FourierTransformer fft;
	double [] phi0, phi0_bar, unstableSoln; // Background stripe configuration and this configuration convoluted with potential.
	double [] sf;
	int ky;
	double kRChunk; //=2piR/L
	public int Lp;

	Grid phiGrid = new Grid("phi field");
	Grid sfGrid = new Grid("SF");
	Plot sfPlot = new Plot("SFs");
	
	Accumulator k0 = new Accumulator();

	public static void main(String[] args) {
		new Control(new Phi4App(), "Svt Ising Field");
	}

	public void load(Control c) {
		c.frameTogether("Grids", phiGrid, sfPlot, sfGrid);
		params.add("Dynamics", new ChoiceValue("Phi4 corberi","Phi4 NCOP","Phi4 COP"));
		params.addm("Zoom", new ChoiceValue("No", "Yes"));
		params.add("Init Conditions", new ChoiceValue("Random Gaussian", "Read 1D Soln","Read From File", "Constant"));
		params.addm("Noise", new DoubleValue(1.0, 0.0, 1.0).withSlider());
		params.addm("T", 0.1);
		params.addm("H", 0.6);
		params.addm("R", 1.0);
		params.addm("Random seed", 0);
		params.addm("dx", 50);
		params.add("Magnetization", 0.0);
		params.addm("Lp", 128);
		params.addm("g", 1.0);
		params.addm("r", -0.000003);
		params.addm("dt", 0.1);
		params.add("mean phi");
		params.add("Time");
		flags.add("Clear");
	}

	public void animate() {
		params.set("Time", field.time());
		params.set("Lp", field.Lp);
		params.set("mean phi", field.mean(field.phi));
		phiGrid.setAutoScale(false);
		phiGrid.registerData(Lp,Lp,field.phi);
		sfGrid.registerData(Lp, Lp, sf);
//		sfPlot.setLogScale(false, true);
		sfPlot.setAutoScale(true);
		sfPlot.registerPoints("sf", k0, Color.BLUE);

		if(flags.contains("Clear")){
			flags.clear();			
		}

	}

	public void clear() {
		k0.clear();
	}

	public void run() {
		clear();
		initialize();
		System.out.println("init");
		String dynamics = params.sget("Dynamics");
		while (true) {
			field.readParams(params);
			if (dynamics == "Phi4 corberi")	
				field.simulatePhi4Corberi();
			else if (dynamics == "Phi4 COP")
				field.simulatePhi4COP();
			else if (dynamics == "Phi4 NCOP")
				field.simulatePhi4NCOP();
			sf = fft.calculate2DSF(field.phi, false, false);
//			k0.accum(field.time(), sf[Lp*Lp/2+Lp/2]);
			k0.accum(field.t, sf[0]);
			Job.animate();

		}	
	}

	private void initialize(){

		field = new phi4(params);
		this.Lp=field.Lp;

		fft = new FourierTransformer(Lp);
		phi0 = new double [Lp];
		unstableSoln = new double [Lp];
	}	


	

}

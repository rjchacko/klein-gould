package rachele.ising.dim2.apps;

import java.awt.Color;
import rachele.ising.dim2.IsingField2D;
import rachele.ising.dim2.StripeClumpFieldSim;
//import rachele.util.FileUtil;
import rachele.util.FourierTransformer;
import scikit.dataset.Accumulator;
import scikit.dataset.PointSet;
import scikit.graphics.dim2.Geom2D;
import scikit.graphics.dim2.Grid;
import scikit.graphics.dim2.Plot;
import scikit.jobs.Control;
import scikit.jobs.Job;
import scikit.jobs.Simulation;
import scikit.jobs.params.ChoiceValue;
import scikit.jobs.params.DirectoryValue;
import scikit.jobs.params.DoubleValue;
import scikit.jobs.params.FileValue;

public class SvtConservedFieldApp extends Simulation{

	IsingField2D ising;
	FourierTransformer fft;
    StripeClumpFieldSim sc;
	double [] eta, etaK, sf; //right hand side
	double [] phi0, phi0_bar, unstableSoln; // Background stripe configuration and this configuration convoluted with potential.
	int ky;
	double kRChunk; //=2piR/L
	public String writeDir;

	//RUN OPTIONS
	boolean writeToFile = false;
	boolean accumsOn = false;
	

	Accumulator sf_hAcc = new Accumulator();
	Accumulator sf_vAcc = new Accumulator();
	
	
	public int Lp;
	Grid phiGrid = new Grid("phi field");
	Grid sfGrid = new Grid("sf");
//	Grid etaDotSF = new Grid("sf phi field");
	Plot plotSF = new Plot("SFs");
	Plot hSlice = new Plot("Horizontal Slice");
	Plot vSlice = new Plot("Vertical Slice"); 


	public static void main(String[] args) {
		new Control(new SvtConservedFieldApp(), "Svt Ising Conserved Field");
	}

	public void load(Control c) {
		c.frameTogether("Grids", phiGrid, sfGrid, hSlice, plotSF);
		params.add("Data Dir",new DirectoryValue("/Users/erdomi/data/lraim/stripeToClumpInvestigation/ftResults/svtFieldApp/stage1/constantH/"));
//		params.add("2D Input File", new FileValue("/home/erdomi/data/lraim/configs/inputConfig"));
		params.add("1D Input File", new FileValue("/Users/erdomi/data/lraim/configs1dAutoConserved/L128R45T0.04m0.7"));
		params.add("New 1D Input File", new FileValue("/Users/erdomi/data/lraim/configs1dAutoConserved/L128R45T0.04m0.7"));
		params.add("Dynamics", new ChoiceValue("Conserved Mob", "Conserved", "Phi4", "Langevin"));
//		params.addm("Zoom", new ChoiceValue("Yes", "No"));
		params.addm("Interaction", new ChoiceValue("Square", "Circle"));
		params.addm("Dynamics?", new ChoiceValue("Langevin No M Convervation", "Langevin Conserve M"));
		params.add("Init Conditions", new ChoiceValue("Random Gaussian","Read 1D Soln", "Read From File", "Constant"));
//		params.addm("Approx", new ChoiceValue( "Modified Dynamics", "None"));
		params.addm("Noise", new DoubleValue(1.0, 0.0, 1.0).withSlider());
		params.addm("Horizontal Slice", new DoubleValue(0.42, 0, 0.9999).withSlider());
		params.addm("Vertical Slice", new DoubleValue(0.5, 0, 0.9999).withSlider());
		params.addm("kR", new DoubleValue(5.135622302, 0.0, 6.0).withSlider());
		params.addm("T", 0.04);
		params.addm("H", 0.0);
//		params.addm("dT", 1.0);
		params.addm("tolerance", 0.01);
		params.addm("J", -1.0);
		params.addm("R", 4600000.0);
		params.addm("Random seed", 0);
		params.add("L/R", 2.782608696);
		params.add("R/dx", 46.0);
		params.add("kR bin-width", 0.1);
		params.add("Magnetization", 0.0);
		params.addm("ky", 2);
		params.addm("kx", 2);
		params.addm("dt", 1.0);
		params.add("mean phi");
		params.add("Time");
		params.add("Lp");
		flags.add("Clear");
	}

	public void animate() {
		params.set("Time", ising.time());
		params.set("Lp", ising.Lp);
		params.set("mean phi", ising.mean(ising.phi));
		phiGrid.setAutoScale(false);
		phiGrid.registerData(Lp,Lp,ising.phi);
		sfGrid.registerData(Lp,Lp,sf);
		hSlice.setAutoScale(true);
		vSlice.setAutoScale(true);
		
		double horizontalSlice = params.fget("Horizontal Slice");
		double verticalSlice = params.fget("Vertical Slice");

		phiGrid.addDrawable(
				Geom2D.line(0, horizontalSlice, 1, horizontalSlice, Color.GREEN));
		phiGrid.addDrawable(
				Geom2D.line(verticalSlice, 0, verticalSlice, 1, Color.BLUE));
		hSlice.registerLines("Slice", ising.getHslice(horizontalSlice), Color.GREEN);
		hSlice.registerLines("phi0", new PointSet(0, 1, phi0) , Color.BLACK);
		hSlice.registerLines("unstable soln", new PointSet(0, 1, unstableSoln) , Color.RED);
		hSlice.registerLines("Slice h", ising.getVslice(verticalSlice), Color.BLUE);
		plotSF.setLogScale(false, true);
//		plotSF.registerPoints("vert", sf_vAcc ,Color.BLUE);
		plotSF.registerPoints("hor", sf_hAcc ,Color.GREEN);
		
		vSlice.registerLines("Slice", ising.getVslice(verticalSlice), Color.BLUE);
		if(flags.contains("Clear")){
			sf_hAcc.clear();
			sf_vAcc.clear();
			flags.clear();			
		}

	}

	public void clear() {
		sf_vAcc.clear();
		sf_hAcc.clear();		
	}

	public void run() {
		clear();
		writeDir = params.sget("Data Dir");
		initialize();
		System.out.println("init");
		for (int i = 0; i < Lp*Lp; i++)
			eta[i] = ising.phi[i] - phi0[i%Lp];
		calcHspinodal();
		String dynamics = params.sget("Dynamics");
		int kyInt = params.iget("ky");
		int kxInt = params.iget("kx");
		
		while (true) {
			ising.readParams(params);
			if (dynamics == "Langevin"){
					ising.simulate();// ising.simulateSimple();
				}else if (dynamics == "Glauber"){
					ising.simulateGlauber();
				}else if (dynamics == "Conserved"){
					ising.simulateConservedFiniteDiff();
				}else if (dynamics == "Phi4"){
					ising.simulatePhi4();
				}else if (dynamics == "Conserved Mob")
					ising.simulateConservedFiniteDiffMob();
			sf = fft.calculate2DSF(ising.phi, false, true);
//			sf_hAcc.accum(ising.time(), sf[kxInt]);
//			sf_hAcc.accum(ising.time(), sf[(Lp - kxInt) % Lp]);
			sf_hAcc.accum(ising.time(),sf[kyInt*Lp + kxInt ]);
			sf_vAcc.accum(ising.time(), sf[kyInt*Lp]);
			sf_vAcc.accum(ising.time(), sf[(Lp - kyInt)*Lp]);
			Job.animate();
		}
	}

	private void calcHspinodal(){
		double rho = Math.sqrt(1+ising.T/(ising.findVkSquare(IsingField2D.KRsquare,0.0)));
		double hs = rho + (ising.T/2.0)*(Math.log(1.0+rho) - Math.log (1-rho));
		System.out.println("H spinodal for T = " + ising.T + " is " + hs);
	}



	private void initialize(){

		ising = new IsingField2D(params);
		this.Lp=ising.Lp;
		sc = new StripeClumpFieldSim(ising, params);
		fft = new FourierTransformer(Lp);
		eta = new double[Lp*Lp];
		etaK = new double[Lp*Lp];
		phi0 = new double [Lp];
		unstableSoln = new double [Lp];
		String fileName = params.sget("1D Input File");
		phi0 = ising.getSymmetricSlice(fileName);
		fileName = params.sget("New 1D Input File");
		unstableSoln = ising.getSymmetricSlice(fileName);
	}	

}	

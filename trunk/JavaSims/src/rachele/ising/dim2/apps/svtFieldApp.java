package rachele.ising.dim2.apps;

//import static scikit.util.Utilities.asList;

import java.awt.Color;
import java.io.File;
import rachele.ising.dim2.IsingField2D;
import rachele.ising.dim2.StripeClumpFieldSim;
import rachele.util.FileUtil;
import rachele.util.FourierTransformer;
import scikit.dataset.Accumulator;
import scikit.dataset.PointSet;
//import scikit.graphics.dim2.Geom2D;
import scikit.graphics.dim2.Grid;
import scikit.graphics.dim2.Plot;
import scikit.jobs.Control;
import scikit.jobs.Job;
import scikit.jobs.Simulation;
import scikit.jobs.params.ChoiceValue;
import scikit.jobs.params.DirectoryValue;
import scikit.jobs.params.DoubleValue;
import scikit.jobs.params.FileValue;

/**
* 
* Something is screwed up with the first run of this program.  
* I don't know why, but you can avoid the problems for some reason
* by starting the program, reseting and starting over.
* Then it seems to be fine.
* 
* This program is just for testing (it has no max time.)  Disable the 
* accumulators to look at long times.
* 
*  Records the structure factors vs time for indefinite time for
*  several k values.
*
*/
public class svtFieldApp extends Simulation{

	IsingField2D ising;
	FourierTransformer fft;
    StripeClumpFieldSim sc;
	double [] eta, etaK, sf; //right hand side
	double [] phi0, phi0_bar; // Background stripe configuration and this configuration convoluted with potential.
	int ky;
	double kRChunk; //=2piR/L
	public String writeDir;

	//RUN OPTIONS
	boolean writeToFile = true;
	boolean accumsOn = true;
	/**
	 * 
	 */
	boolean calcContribs = true;
	int accNo = 1;
	
	Accumulator [] etaAcc = new Accumulator [accNo];
	Accumulator [] etaLTAcc = new Accumulator [accNo];
	Accumulator [] etaLTkAcc = new Accumulator [accNo];

    int [][][] sfLabel = new int [2][4][accNo];
	
	public int Lp;
	Grid etaDot = new Grid("ising delta phi");
	Grid phiGrid = new Grid("phi field");
//	Grid etaDotSF = new Grid("sf phi field");
	Plot SFvTime = new Plot("svt");
	Plot EtavTime = new Plot("etavt");
	Plot eta2 = new Plot("etavt2");
//	Grid phi0Grid = new Grid("phi_0");
//	Grid phi0SliceGrid = new Grid("phi_0 slice");
//	Plot etaVsTimeLC = new Plot("eta v t LC");
//	Plot fkPlot = new Plot ("f(k)");
//	Plot etakSimPlot = new Plot("eta k sim");
//	Plot convolve = new Plot("convolve");
	Plot hSlice = new Plot("Horizontal Slice");
	Plot vSlice = new Plot("Vertical Slice"); 
	Grid mobilityCheck = new Grid("Mobility");
//	Plot eVector = new Plot ("Largest Eigenvector");


	public static void main(String[] args) {
		new Control(new svtFieldApp(), "Svt Ising Field");
	}

	public void load(Control c) {
		c.frameTogether("Grids", phiGrid, vSlice, hSlice, mobilityCheck);
		c.frame(EtavTime);
		c.frameTogether("etas", etaDot, eta2);
		params.add("Data Dir",new DirectoryValue("/home/erdomi/data/lraim/stripeToClumpInvestigation/ftResults/svtFieldApp/testRuns"));
		params.add("2D Input File", new FileValue("/home/erdomi/data/lraim/configs/inputConfig"));
		params.add("1D Input File", new FileValue("/home/erdomi/data/lraim/configs1dAutoName/L128R45T0.04h0.8"));
		params.add("New 1D Input File", new FileValue("/home/erdomi/data/lraim/configs1dAutoName/L128R45T0.04h0.8"));
		params.add("Dynamics", new ChoiceValue("Langevin", "Glauber" ));
		params.addm("Zoom", new ChoiceValue("Yes", "No"));
		params.addm("Interaction", new ChoiceValue("Square", "Circle"));
		params.addm("Dynamics?", new ChoiceValue("Langevin No M Convervation", "Langevin Conserve M"));
		params.add("Init Conditions", new ChoiceValue("Read 1D Soln", "Read From File","Random Gaussian", "Constant"));
		params.addm("Approx", new ChoiceValue("None", "Modified Dynamics"));
		params.addm("Noise", new DoubleValue(0.0, 0.0, 1.0).withSlider());
		params.addm("Horizontal Slice", new DoubleValue(0.5, 0, 0.9999).withSlider());
		params.addm("Vertical Slice", new DoubleValue(0.5, 0, 0.9999).withSlider());
		params.addm("kR", new DoubleValue(5.135622302, 0.0, 6.0).withSlider());
		params.addm("T", 0.04);
		params.addm("H", 0.80);
		params.addm("dT", 0.001);
		params.addm("tolerance", 0.01);
		params.addm("J", -1.0);
		params.addm("R", 2000000.0);
		params.addm("Random seed", 0);
		params.add("L/R", 2.782608696);
		params.add("R/dx", 50.0);
		params.add("kR bin-width", 0.1);
		params.add("Magnetization", 0.0);
		params.addm("ky", 2);
		params.addm("dkx", 1);
		params.addm("dt", 0.01);
		params.add("mean phi");
		params.add("Time");
		params.add("Lp");
		flags.add("Clear");
		flags.add("Write 1D Config");
	}

	public void animate() {
		params.set("Time", ising.time());
		params.set("Lp", ising.Lp);
		params.set("mean phi", ising.mean(ising.phi));
		phiGrid.registerData(Lp, Lp,(fft.calculate2DSF(ising.phi, true, true)));
		etaDot.setAutoScale(true);
		etaDot.registerData(Lp, Lp, ising.phi);
//		System.out.println("thel " + sc.eta_bar2[Lp]);
//		etaDotSF.registerData(Lp, Lp, etaK);
//		etaVsTimeLC.setAutoScale(true);
		hSlice.setAutoScale(true);
		vSlice.setAutoScale(true);
//		eVector.setAutoScale(true);
		SFvTime.setAutoScale(true);
		SFvTime.setLogScale(false, true);
		EtavTime.setAutoScale(true);
		EtavTime.setLogScale(false, true);
		eta2.setAutoScale(true);
		eta2.setLogScale(false, true);
		
		double horizontalSlice = params.fget("Horizontal Slice");
		double verticalSlice = params.fget("Vertical Slice");

//		phiGrid.setDrawables(asList(
//				Geom2D.line(0, horizontalSlice, 1, horizontalSlice, Color.GREEN),
//				Geom2D.line(verticalSlice, 0, verticalSlice, 1, Color.BLUE)));

		mobilityCheck.registerData(Lp, Lp, sc.etaBar_k);
		
		hSlice.registerLines("Slice", ising.getHslice(horizontalSlice), Color.GREEN);
		hSlice.registerLines("phi0", new PointSet(0, 1, phi0) , Color.BLACK);
		vSlice.registerLines("Slice", ising.getVslice(verticalSlice), Color.BLUE);
		if(accumsOn){
			for (int i = 0; i < accNo; i ++){
//				float colorChunk = (float)i/(float)accNo;
//				Color col = Color.getHSBColor(colorChunk, 1.0f, 1.0f);
				StringBuffer sb = new StringBuffer();sb.append("s(t) Ave "); sb.append(i);
				EtavTime.registerLines(sb.toString(), etaAcc[i], Color.black);
				StringBuffer sb2 = new StringBuffer();sb2.append("etaLT "); sb2.append(i);
				EtavTime.registerLines(sb2.toString(), etaLTAcc[i], Color.BLUE);

				StringBuffer sb3 = new StringBuffer();sb3.append("etaLT_k "); sb3.append(i);
				
				EtavTime.registerLines(sb3.toString(), etaLTkAcc[i], Color.RED);
				eta2.registerLines(sb3.toString(), etaLTkAcc[i], Color.RED);
				
//				etaDot.registerData(Lp, Lp, sc.etaKchange);
			}
		}
		if(flags.contains("Write 1D Config"))
			write1Dconfig();
		if(flags.contains("Clear")){
			flags.clear();			
		}

	}

	public void clear() {
	}

	public void run() {
		clear();
		writeDir = params.sget("Data Dir");
		initialize();
		System.out.println("init");
		String approx = params.sget("Approx");
		double recordStep = 0.00001;	
		for (int i = 0; i < Lp*Lp; i++)
			eta[i] = ising.phi[i] - phi0[i%Lp];
		sc.initEta();
		calcHspinodal();
		Job.animate();
		initFiles();
		String dynamics = params.sget("Dynamics");

		
		while (true) {
			ising.readParams(params);
			if(calcContribs){ 
				if(dynamics == "Langevin"){
					if(approx == "None"){
						ising.simCalcContrib();
						recordContribToFile();
						sc.simulateLinear();
						sc.simulateLinearK();
//						sc.simulateLinearKbar1d();
					}else if(approx == "Modified Dynamics"){
						ising.simModCalcContrib();
						recordContribToFile();	
					}
				}else if(dynamics == "Glauber"){
					ising.glauberCalcContrib();
					recordContribToFile();
//					System.out.println("Gl calc cont");
				}
			}
			else{
				if (dynamics == "Langevin"){
					if (approx == "None") ising.simulate();// ising.simulateSimple();
					else if (approx == "Modified Dynamics")ising.simulate();
				}else if (dynamics == "Glauber")
					ising.simulateGlauber();
			}
			if(ising.time() >= recordStep){
				for (int i = 0; i < Lp*Lp; i++)
					eta[i] = ising.phi[i] - phi0[i%Lp];
				etaK = fft.calculate2DFT(eta);
				sf = fft.find2DSF(ising.phi, ising.L);
				double [] scEtaK =fft.calculate2DFT(sc.etaLT); 
//				double [] scEtaK2 =fft.find2DSF(sc.etaLT2,ising.L); 
				if(accumsOn){
					for (int i = 0; i < accNo; i++){
						etaAcc[i].accum(ising.time(), etaK[sfLabel[0][1][i]]);
						etaLTAcc[i].accum(ising.time(), scEtaK[sfLabel[0][1][i]]);	
						int dkx = params.iget("dkx");
//						etaLTkAcc[i].accum(ising.time(), scEtaK2[sfLabel[0][1][i]]);
//						etaLTkAcc[i].accum(ising.time(), Math.pow(sc.etaLT2D_k[sfLabel[0][1][i]],1));
						etaLTkAcc[i].accum(ising.time(), Math.pow(sc.etaLT_k[i*dkx],1));
//						System.out.println(Math.pow(sc.etaLT_k[i*dkx],1));
					}
				}
				recordSfDataToFile(sf);
				recordStep += 0.001;
			}

			Job.animate();
//			System.out.println("conserve " + ising.magConservation);

		}	
	}

	private void calcHspinodal(){
		double rho = Math.sqrt(1+ising.T/(ising.findVkSquare(IsingField2D.KRsquare,0.0)));
		double hs = rho + (ising.T/2.0)*(Math.log(1.0+rho) - Math.log (1-rho));
		System.out.println("H spinodal for T = " + ising.T + " is " + hs);
	}

	private void recordSfDataToFile(double [] data1){		
		String fileName = params.sget("Data Dir") + File.separator + "e0";
		StringBuffer fileBuffer = new StringBuffer(); fileBuffer.append(fileName);
		for (int i=0; i < accNo; i ++){	
			fileBuffer.deleteCharAt(fileBuffer.length()-1);
			fileBuffer.deleteCharAt(fileBuffer.length()-1);
			fileBuffer.append("h");fileBuffer.append(i); fileName = fileBuffer.toString();
			double sum = 0;
			for (int j=0; j<4; j++) sum += data1[sfLabel[0][j][i]]; 
			sum /=4;
			FileUtil.printlnToFile(fileName, ising.time(), sum);	
			fileBuffer.deleteCharAt(fileBuffer.length()-1);
			fileBuffer.deleteCharAt(fileBuffer.length()-1);
			fileBuffer.append("v");fileBuffer.append(i); fileName = fileBuffer.toString();
			sum = 0;
			for (int j=0; j<4; j++) sum += data1[sfLabel[1][j][i]]; 
			sum /=4;
			FileUtil.printlnToFile(fileName, ising.time(), sum);			
		}
	}	

	private void recordContribToFile(){
		String fileName = params.sget("Data Dir") + File.separator + "n";
		double driftC = ising.driftContrib;
		double noiseC = 1.0-driftC;
		FileUtil.printlnToFile(fileName, ising.time(), noiseC);
		fileName = params.sget("Data Dir") + File.separator + "d";
		FileUtil.printlnToFile(fileName, ising.time(), driftC);
	}
	
	private void initFiles(){
		String message1 = "#Field theory of SF vs time for several values of k.  Stripe to clump. ";
		String fileName = params.sget("Data Dir") + File.separator + "e0";
		StringBuffer fileBuffer = new StringBuffer(); fileBuffer.append(fileName);
		for (int i=0; i < accNo; i ++){
			StringBuffer mb = new StringBuffer();
			int label = params.iget("dkx")*i;
			mb.append("#kR line value is ");
			double kRLine = 2*ising.R*Math.PI*(params.iget("ky"))/ising.L; mb.append(kRLine);
			mb.append(" sf label is ");	mb.append(label); mb.append(" kR value = ");
			double krvalue = 2*ising.R*Math.PI*(label)/ising.L;
			mb.append(krvalue);	
			fileBuffer.deleteCharAt(fileBuffer.length()-1);	fileBuffer.deleteCharAt(fileBuffer.length()-1);
			fileBuffer.append("h");fileBuffer.append(i); fileName = fileBuffer.toString();
			String message2 = mb.toString();
			FileUtil.initFile(fileName, params, message1, message2);
			fileBuffer.deleteCharAt(fileBuffer.length()-1);	fileBuffer.deleteCharAt(fileBuffer.length()-1);
			fileBuffer.append("v");fileBuffer.append(i); fileName = fileBuffer.toString();
			FileUtil.initFile(fileName, params, message1, message2);

		}
		if(calcContribs){
			String mes1 = "# Noise contribution to dynamics";
			String noiseFile =  params.sget("Data Dir") + File.separator + "n";
			FileUtil.initFile(noiseFile, params, mes1);			
			String mes2 = "# Drift contribution to dynamics";
			String driftFile =  params.sget("Data Dir") + File.separator + "d";
			FileUtil.initFile(driftFile, params, mes2);			
		}
	}

	private void initialize(){

		ising = new IsingField2D(params);
		this.Lp=ising.Lp;
		sc = new StripeClumpFieldSim(ising, params);
		for (int i = 0; i < accNo; i++){ 
			if(accumsOn){
				etaAcc[i] = new Accumulator();
				etaLTAcc[i] = new Accumulator();
				etaLTkAcc[i] = new Accumulator();
			}
			int index, x, y;
			int dkx = params.iget("dkx");
			ky = params.iget("ky");
			int kx = i*dkx;
			//Horizontal labels
			int size = Lp*Lp;
			y = ky; x = kx;//Up, right
			index = Lp * y+x; 
			sfLabel[0][0][i] = index;
//			System.out.println("kx= "+kx+" ky= "+ky+" Lp = "+Lp+" index = "+index);
			x = (Lp - kx) % Lp;//Up, left	
			index = Lp * y+x; sfLabel[0][1][i] = index%size;
			y = (Lp - ky) % Lp; x = kx;//Down, right
			index = Lp*y+x; sfLabel[0][2][i] = index%size;
			x = (Lp - kx) % Lp;//Up, left	
			index = Lp * y + x; sfLabel[0][3][i] = index%size;
			//Vertical labels
			x = ky; y = kx;//Right, Up
			index = Lp * y + x; sfLabel[1][0][i] = index%size;
			y = (Lp - kx)%Lp;//Right, Down	
			index = Lp * y + x; sfLabel[1][1][i] = index%size;
			x = (Lp - ky ) % Lp; y = kx;//Left, Up
			index = Lp * y + x; sfLabel[1][2][i] = index%size;
			y = (Lp - kx) % Lp;//Left, Down
			index = Lp * y + x; sfLabel[1][3][i] = index%size;

		}

		fft = new FourierTransformer(Lp);
		eta = new double[Lp*Lp];
		etaK = new double[Lp*Lp];
		phi0 = new double [Lp];
		String fileName = params.sget("1D Input File");
		phi0 = ising.getSymmetricSlice(fileName);
	}	

    private void write1Dconfig(){
        String configFileName = params.sget("New 1D Input File");
        FileUtil.deleteFile(configFileName);
        FileUtil.writeConfigToFile(configFileName, ising.Lp, ising.phi);
        System.out.println("Config written");
        flags.clear();
}

	
}
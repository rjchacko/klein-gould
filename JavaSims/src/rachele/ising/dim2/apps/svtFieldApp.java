package rachele.ising.dim2.apps;

import static scikit.util.Utilities.asList;

import java.awt.Color;
import java.io.DataInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import rachele.ising.dim2.IsingField2D;
import rachele.util.FileUtil;
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
public class svtFieldApp extends Simulation{

	IsingField2D ising;
	FourierTransformer fft;
	double [] eta, etaK, sf; //right hand side
	double [] phi0, phi0_bar; // Background stripe configuration and this configuration convoluted with potential.
	int ky;
	double kRChunk; //=2piR/L
	boolean clearFile;
	public String writeDir;

	//RUN OPTIONS
	boolean writeToFile = true;
	
	int accNo = 6;
    Accumulator [] etaAcc = new Accumulator [accNo];
    Accumulator [] sfAcc = new Accumulator [accNo];
    int [] sfLabel = new int [accNo];
	
	public int Lp;
	Grid etaDot = new Grid("ising delta phi");
	Grid phiGrid = new Grid("phi field");
	Grid etaDotSF = new Grid("sf phi field");
	Plot SFvTime = new Plot("svt");
	Plot EtavTime = new Plot("etavt");
	Grid phi0Grid = new Grid("phi_0");
	Grid phi0SliceGrid = new Grid("phi_0 slice");
	Plot etaVsTimeLC = new Plot("eta v t LC");
	Plot fkPlot = new Plot ("f(k)");
	Plot etakSimPlot = new Plot("eta k sim");
	Plot convolve = new Plot("convolve");
	Plot hSlice = new Plot("Horizontal Slice");
	Plot vSlice = new Plot("Vertical Slice"); 
	Plot eVector = new Plot ("Largest Eigenvector");


	public static void main(String[] args) {
		new Control(new svtFieldApp(), "Ising Linear Test");
	}

	public void load(Control c) {
		c.frameTogether("Grids", phiGrid, vSlice, hSlice);
		params.add("Data Dir",new DirectoryValue("/home/erdomi/data/lraim/stripeToClumpInvestigation/ftResults/svtFieldApp"));
		params.add("2D Input File", new FileValue("/home/erdomi/data/lraim/configs/inputConfig"));
		params.add("1D Input File", new FileValue("/home/erdomi/data/lraim/configs1d/L128R46T0-04h0-8"));
		params.add("1D phi0 File", new FileValue("/home/erdomi/data/lraim/configs1d/L128R46T0-04h0-8"));
		params.addm("Zoom", new ChoiceValue("Yes", "No"));
		params.addm("Interaction", new ChoiceValue("Square", "Circle"));
		params.addm("Dynamics?", new ChoiceValue("Langevin No M Convervation"));
		params.add("Init Conditions", new ChoiceValue("Read 1D Soln", "Read From File","Random Gaussian"));
		params.addm("Approx", new ChoiceValue("Slow", "HalfStep", "TimeAdjust", "Phi4","Phi4HalfStep"));
		params.addm("Noise", new DoubleValue(1.0, 0.0, 1.0).withSlider());
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
		params.add("L/R", 2.7826087);
		params.add("R/dx", 50.0);
		params.add("kR bin-width", 0.1);
		params.add("Magnetization", 0.0);
		params.addm("ky", 2);
		params.addm("dt", 0.001);
		params.add("Time");
		//params.add("Mean Phi");
		params.add("Lp");
		flags.add("Clear");
		flags.add("Write 1D Config");
	}

	public void animate() {
		params.set("Time", ising.time());
		params.set("Lp", ising.Lp);
		phiGrid.registerData(Lp, Lp, ising.phi);
		etaDotSF.registerData(Lp, Lp, etaK);
		etaVsTimeLC.setAutoScale(true);
		hSlice.setAutoScale(true);
		vSlice.setAutoScale(true);
		eVector.setAutoScale(true);
		SFvTime.setAutoScale(true);
		SFvTime.setLogScale(false, true);
		EtavTime.setAutoScale(true);
		EtavTime.setLogScale(false, true);
		
		double horizontalSlice = params.fget("Horizontal Slice");
		double verticalSlice = params.fget("Vertical Slice");

		phiGrid.setDrawables(asList(
				Geom2D.line(0, horizontalSlice, 1, horizontalSlice, Color.GREEN),
				Geom2D.line(verticalSlice, 0, verticalSlice, 1, Color.BLUE)));

		hSlice.registerLines("Slice", ising.getHslice(horizontalSlice), Color.GREEN);
		hSlice.registerLines("phi0", new PointSet(0, 1, phi0) , Color.BLACK);
		vSlice.registerLines("Slice", ising.getVslice(verticalSlice), Color.BLUE);
		for (int i = 0; i < accNo; i ++){
			float colorChunk = (float)i/(float)accNo;
			Color col = Color.getHSBColor(colorChunk, 1.0f, 1.0f);
			StringBuffer sb = new StringBuffer();sb.append("s(t) Ave "); sb.append(i);
			EtavTime.registerLines(sb.toString(), etaAcc[i], col);
			StringBuffer sb2 = new StringBuffer();sb2.append("s(t) "); sb2.append(i);
			SFvTime.registerLines(sb2.toString(), sfAcc[i], col);
		}
		
		if(flags.contains("Clear")){
			flags.clear();			
		}

	}

	public void clear() {
	}

	public void run() {
		writeDir = params.sget("Data Dir");
		if (flags.contains("Write 1D Config")){
			write1Dconfig();
			flags.clear();
		}
		clearFile = true;
		initialize();
		System.out.println("init");
		ky = params.iget("ky");

		double recordStep = 0.00001;	

		for (int i = 0; i < Lp*Lp; i++)
			eta[i] = ising.phi[i] - phi0[i%Lp];

		System.out.println("ky = " + ky);
		calcHspinodal();
		Job.animate();
		initFiles();
				
		while (true) {
			ising.readParams(params);
			//ising.simulate();
			ising.simulateSimple();
			if(ising.time() >= recordStep){
				for (int i = 0; i < Lp*Lp; i++)
					eta[i] = ising.phi[i] - phi0[i%Lp];
				etaK = fft.find2DSF(eta, ising.L);
				sf = fft.find2DSF(ising.phi, ising.L);
				for (int i = 0; i < accNo; i++){
					etaAcc[i].accum(ising.time(), etaK[sfLabel[i]]);
					sfAcc[i].accum(ising.time(), sf[sfLabel[i]]);					
				}
				recordSfDataToFile(etaK, sf);
				recordStep += .0001;
			}

			Job.animate();

		}	
	}

	private void calcHspinodal(){
		double rho = Math.sqrt(1+ising.T/(ising.findVkSquare(IsingField2D.KRsquare,0.0)));
		double hs = rho + (ising.T/2.0)*(Math.log(1.0+rho) - Math.log (1-rho));
		System.out.println("H spinodal for T = " + ising.T + " is " + hs);
	}

	public void recordSfDataToFile(double [] data1, double [] data2){
		
		String fileName = params.sget("Data Dir") + File.separator + "e0";
		StringBuffer fileBuffer = new StringBuffer(); fileBuffer.append(fileName);
		for (int i=0; i < accNo; i ++){
			
			StringBuffer mb = new StringBuffer();
			mb.append("# sf label = ");	mb.append(sfLabel[i]); mb.append(" kR value = ");
			double krvalue = 2*ising.R*Math.PI*(sfLabel[i])/ising.L;
			mb.append(krvalue);			
			fileBuffer.deleteCharAt(fileBuffer.length()-1);	fileBuffer.deleteCharAt(fileBuffer.length()-1);
			fileBuffer.append("e");fileBuffer.append(i); fileName = fileBuffer.toString();
			FileUtil.printlnToFile(fileName, ising.time(), data1[sfLabel[i]]);	
			fileBuffer.deleteCharAt(fileBuffer.length()-1);	fileBuffer.deleteCharAt(fileBuffer.length()-1);
			fileBuffer.append("s");fileBuffer.append(i); fileName = fileBuffer.toString();
			FileUtil.printlnToFile(fileName, ising.time(), data2[sfLabel[i]]);			
		}
	}	

	private void initFiles(){
		String message1 = "#Field theory of SF vs time for several values of k.  Stripe to clump. ";
		String fileName = params.sget("Data Dir") + File.separator + "e0";
		StringBuffer fileBuffer = new StringBuffer(); fileBuffer.append(fileName);
		for (int i=0; i < accNo; i ++){
			StringBuffer mb = new StringBuffer();
			mb.append("# sf label = ");	mb.append(sfLabel[i]); mb.append(" kR value = ");
			double krvalue = 2*ising.R*Math.PI*(sfLabel[i])/ising.L;
			mb.append(krvalue);			
			fileBuffer.deleteCharAt(fileBuffer.length()-1);	fileBuffer.deleteCharAt(fileBuffer.length()-1);
			fileBuffer.append("e");fileBuffer.append(i); fileName = fileBuffer.toString();
			String message2 = mb.toString();
			FileUtil.initFile(fileName, params, message1, message2);
			fileBuffer.deleteCharAt(fileBuffer.length()-1);	fileBuffer.deleteCharAt(fileBuffer.length()-1);
			fileBuffer.append("s");fileBuffer.append(i); fileName = fileBuffer.toString();
			FileUtil.initFile(fileName, params, message1, message2);		
		}
	}
	
	public void readInputParams(String FileName){
		try {
			File inputFile = new File(FileName);
			DataInputStream dis = new DataInputStream(new FileInputStream(inputFile));
			double readData;

			readData = dis.readDouble();
			System.out.println(readData);
			params.set("J", readData);
			dis.readChar();				

			readData = dis.readDouble();
			System.out.println(readData);
			params.set("R", readData);
			dis.readChar();	

			readData = dis.readDouble();
			System.out.println(readData);
			params.set("L/R", readData);
			dis.readChar();

			readData = dis.readDouble();
			System.out.println(readData);
			params.set("R/dx", readData);
			dis.readChar();	

			System.out.println("input read");
		}catch(IOException ex){
			ex.printStackTrace();
		}
	}

	void initialize(){
		if(params.sget("Init Conditions") == "Read From File")
			readInputParams(params.sget("2D Input File"));
		ising = new IsingField2D(params);
		for (int i = 0; i < accNo; i++){
			etaAcc[i] = new Accumulator();
			sfAcc[i] = new Accumulator();
			sfLabel[i] = ky*Lp+i;
		}
		this.Lp = ising.Lp;
		fft = new FourierTransformer(Lp);
		eta = new double[Lp*Lp];
		etaK = new double[Lp*Lp];
		phi0 = new double [Lp];
		String fileName = params.sget("1D phi0 File");
		phi0 = ising.getSymmetricSlice(fileName);
	}	

	private void write1Dconfig(){
		String configFileName = params.sget("1D Input File");
		FileUtil.deleteFile(configFileName);
		FileUtil.writeConfigToFile(configFileName, ising.Lp, ising.phi);
	}
}
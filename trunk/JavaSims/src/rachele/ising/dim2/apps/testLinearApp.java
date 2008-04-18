package rachele.ising.dim2.apps;

import static java.lang.Math.floor;
import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.awt.Color;
import rachele.ising.dim2.FindCoefficients;
import rachele.ising.dim2.FourierTransformer;
import rachele.ising.dim2.IsingField2D;
import rachele.util.FileUtil;
import scikit.graphics.dim2.Grid;
import scikit.graphics.dim2.Plot;
import scikit.jobs.Control;
import scikit.jobs.Job;
import scikit.jobs.Simulation;
import scikit.jobs.params.ChoiceValue;
import scikit.jobs.params.DoubleValue;
//import scikit.numerics.fft.managed.ComplexDouble2DFFT;
import scikit.dataset.PointSet;

public class testLinearApp extends Simulation{
    IsingField2D ising;
    public FindCoefficients coeff;
    public FourierTransformer fft;
    public double [] eta_dot, eta_dot_x;
    public double [] phiFT, del_phi_k, linearConfig;
    public double [] rhs;  // right hand side of equation
    public double kRChunk; //=2piR/L
    Grid etaDotLeftGrid = new Grid("ising delta phi");
    Grid etaDotRightGrid = new Grid("linear theory");
    Grid phiGrid = new Grid("phi field");
    Grid ftPhiField = new Grid("sf phi field");
    Grid phiDotGrid = new Grid("sf linear theory");
    Grid delPhiKGrid = new Grid("delPhi_k");
    Plot slicePlot1 = new Plot("slices");
    Plot slicePlot2 = new Plot("slices");
    
	public static void main(String[] args) {
		new Control(new testLinearApp(), "Ising Linear Test");
	}
	public void load(Control c) {
	c.frameTogether("Grids", etaDotLeftGrid, delPhiKGrid,phiGrid,etaDotRightGrid,ftPhiField  );
	//c.frameTogether("slices", slicePlot1, slicePlot2);
	slicePlot1.setAutoScale(true);
	slicePlot2.setAutoScale(true);
	
	params.addm("Zoom", new ChoiceValue("Yes", "No"));
	params.addm("Interaction", new ChoiceValue("Square", "Circle"));
	params.addm("Dynamics?", new ChoiceValue("Langevin No M Convervation"));
	params.add("Init Conditions", new ChoiceValue("Read 1D Soln", "Read From File","Random Gaussian"));
	params.addm("Approx", new ChoiceValue("Slow", "HalfStep", "TimeAdjust", "Phi4","Phi4HalfStep"));
	params.addm("Noise", new DoubleValue(0.0, 0.0, 1.0).withSlider());
	params.addm("Horizontal Slice", new DoubleValue(0.5, 0, 0.9999).withSlider());
	params.addm("Vertical Slice", new DoubleValue(0.5, 0, 0.9999).withSlider());
	params.addm("kR", new DoubleValue(5.135622302, 0.0, 6.0).withSlider());
	params.addm("T", 0.04);
	params.addm("H", 0.80);
	params.addm("dT", 0.001);
	params.addm("tolerance", 0.0001);
	params.addm("dt", 1.0);
	params.addm("J", -1.0);
	params.addm("R", 2000000.0);
	params.addm("Random seed", 0);
	params.add("L/R", 3.0);
	params.add("R/dx", 43.0);
	params.add("kR bin-width", 0.1);
	params.add("Magnetization", 0.0);
	params.add("Time");
	params.add("Mean Phi");
	params.add("Lp");
	}

	public void animate() {
		params.set("Time", ising.time());
		params.set("Mean Phi", ising.mean(ising.phi));
		params.set("Lp", ising.Lp);
//		for(int i = 0; i < ising.Lp; i ++){
//			eta_dot[i]=0;
//			del_phi_k[i] = 0;
//		}
		etaDotLeftGrid.registerData(ising.Lp, ising.Lp, ising.delPhi);
		delPhiKGrid.registerData(ising.Lp, ising.Lp, del_phi_k);
		etaDotRightGrid.registerData(ising.Lp, ising.Lp, rhs);
		phiGrid.registerData(ising.Lp, ising.Lp, ising.phi);
		ftPhiField.registerData(ising.Lp, ising.Lp, phiFT);
		phiDotGrid.registerData(ising.Lp, ising.Lp, ising.delPhi);
		double [] etaDotSlice = new double [ising.Lp];
		double [] rhsSlice = new double [ising.Lp];
		for (int i = 0; i < ising.Lp; i ++){
			etaDotSlice[i] = eta_dot[i];
			rhsSlice[i] = rhs[i];
		}
		slicePlot1.registerLines("eta dot slice", new PointSet(0, 1, etaDotSlice), Color.BLUE);
		slicePlot2.registerLines("eta dot slice", new PointSet(0, 1, etaDotSlice), Color.BLUE);
		slicePlot2.registerLines("rhs slice", new PointSet(0, 1, rhsSlice), Color.RED);
	}

	public void clear() {
		
	}

	public void run() {
		if(params.sget("Init Conditions") == "Read From File")
			readInputParams("../../../research/javaData/configs/inputParams");
		ising = new IsingField2D(params);
		coeff = new FindCoefficients(ising.Lp);
		fft = new FourierTransformer(ising.Lp);
		double binWidth = params.fget("kR bin-width");
		binWidth = IsingField2D.KR_SP / floor(IsingField2D.KR_SP/binWidth);
		eta_dot = new double[ising.Lp*ising.Lp];
		eta_dot_x = new double[ising.Lp*ising.Lp];
		phiFT = new double[ising.Lp*ising.Lp];
		rhs = new double[ising.Lp*ising.Lp];
		del_phi_k = new double[ising.Lp*ising.Lp];
		linearConfig = new double[ising.Lp*ising.Lp];
		double [] phi_k_old = new double[ising.Lp*ising.Lp];
		double [] phi_k_new = new double[ising.Lp*ising.Lp];
		double [] previousConfig = new double[ising.Lp*ising.Lp];
		double [] eta_x = new double[ising.Lp*ising.Lp];

		for (int i = 0; i < ising.Lp*ising.Lp; i ++)
			previousConfig[i] = ising.phi[i];
		ising.simulate();
		kRChunk = 2*Math.PI/params.fget("L/R");
		
        while (true) {
        	double [] background = new double [ising.Lp*ising.Lp];
        	double [] backgroundSlice = new double [ising.Lp];
        	for (int i = 0; i < ising.Lp*ising.Lp; i++){
        		if(i<ising.Lp)
        			backgroundSlice[i] = previousConfig[i%ising.Lp];
        		background[i] = previousConfig[i];
        	}
        	
    		phi_k_old = fft.calculate2DFT(ising.phi);
    		for (int i = 0; i < ising.Lp*ising.Lp; i ++)
    			eta_x[i] = ising.phi[i]-background[i];
    		//findRightSide(eta_x, previousSlice);
    		//findRightSideNoAssumptions(eta_x, previousConfig);
    		//String fileName = "../../../research/javaData/configs1d/config";
    		//rhs = ising.simulateLinear(eta_x, FileUtil.readConfigFromFile(fileName, ising.Lp));
    		//rhs = ising.simulateLinear(eta_x, background);
    		rhs = ising.simulateLinearWithModDynamics();
    		for (int i = 0; i < ising.Lp*ising.Lp; i ++)
    			previousConfig[i] = ising.phi[i];
    		ising.simulate();
    		phi_k_new = fft.calculate2DFT(ising.phi);
    		for (int i = 0; i < ising.Lp*ising.Lp; i ++){
				eta_dot[i] = (phi_k_old[i] - phi_k_new[i])/ising.dt;
				phiFT[i] = phi_k_new[i];
				eta_dot_x[i] = ising.phi[i]-background[i] -eta_x[i];
    		}
    		del_phi_k = fft.calculate2DFT(ising.delPhi);
			Job.animate();
		}
		
	}
	

	public void findRightSideNoAssumptions(double [] eta_x, double [] config){
		double [] eta_k = new double [ising.Lp*ising.Lp];
		eta_k = fft.calculate2DFT(eta_x);
		coeff.findCoefficientsFromConfig(config);
		for(int i = 0; i < ising.Lp*ising.Lp; i++){
			int currentKx = i%ising.Lp;
			int currentKy = i/ising.Lp;
			double kRxValue = kRChunk*currentKx;
			double kRyValue = kRChunk*currentKy;
			double V_k = findVk(kRxValue, kRyValue);
			double sum = -V_k*eta_k[i];
			for (int j = 0; j < ising.Lp*ising.Lp; j ++){
				int shiftedKx = j%ising.Lp;
				int shiftedKy = j/ising.Lp;
				double coefficient = coeff.aa[shiftedKx][shiftedKy];
				int etaXindex = (currentKx+shiftedKx)%ising.Lp;
				int etaYindex = (currentKy+shiftedKy)%ising.Lp;
				int	etaIndex = etaYindex*ising.Lp + etaXindex;
				sum -= ising.T*coefficient*eta_k[etaIndex];
			}
			rhs[i] = sum;
		}
	}
	
	public void findRightSide(double [] eta_x, double [] slice){
		double [] eta_k = new double [ising.Lp*ising.Lp];
		eta_k = fft.calculate2DFT(eta_x);
		//coeff.findCoefficientsFromSlice(slice);
		coeff.findCoefficientsFromFile();
		for(int i = 0; i < ising.Lp*ising.Lp; i++){
			//coeff.findCoefficientsFromSlice(slice);
			//coeff.findCoefficientsFromFile();
			int kx = i%ising.Lp;
			int ky = i/ising.Lp;
			double kRxValue = kRChunk*kx;
			double kRyValue = kRChunk*ky;
			double V_k = findVk(kRxValue, kRyValue);
			double sum = -V_k*eta_k[i];
			for (int j = 0; j + kx < ising.Lp; j++)
				sum -= ising.T*coeff.a[j]*eta_k[ i + j ];
			for (int j = ising.Lp - kx; j < ising.Lp; j++)
				sum -= ising.T*coeff.a[j]*eta_k[ i + j - ising.Lp ];
			rhs[i] = sum;
		}
	}
	
	public double findVk(double kRx, double kRy){
		double Vkx = (kRx == 0) ? 1 : Math.sin(kRx)/kRx;
		double Vky = (kRy == 0) ? 1 : Math.sin(kRy)/kRy;
		double Vk = Vkx*Vky;
		return Vk;
	}
	
//	public double [] calculate2DFT(double [] data){
//		ComplexDouble2DFFT fft = new ComplexDouble2DFFT(ising.Lp, ising.Lp);
//		double [] fftScratch = new double [2*ising.Lp*ising.Lp];
//		for (int i = 0; i < ising.Lp*ising.Lp; i ++){
//			fftScratch[2*i] = data[i];
//			fftScratch[2*i+1] = 0;
//		}
//		fft.transform(fftScratch);
//		double [] dataFT = new double [ising.Lp*ising.Lp];
//		for (int i = 0; i < ising.Lp*ising.Lp; i ++)
//			dataFT[i] = fftScratch[2*i];
//		return dataFT;
//	}
	
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
	
	public void writeConfiguration(){
		String configFileName = "../../../research/javaData/configs/inputConfig";
		String inputFileName = "../../../research/javaData/configs/inputParams";
		FileUtil.deleteFile(configFileName);
		FileUtil.deleteFile(inputFileName);
		writeInputParams(inputFileName);	
		writeConfigToFile(configFileName);
		System.out.println("Config writtern to file");
	}
	
	public void writeInputParams(String FileName){
		try {
			File inputFile = new File(FileName);
			DataOutputStream dos = new DataOutputStream(new FileOutputStream(inputFile, true));
			
			dos.writeDouble(params.fget("J"));
			dos.writeChar('\t');
			dos.writeDouble(params.fget("R"));
			dos.writeChar('\t');
			dos.writeDouble(params.fget("L/R"));
			dos.writeChar('\t');
			dos.writeDouble(params.fget("R/dx"));
			dos.writeChar('\t');
			dos.writeChar('\n');
			dos.close();
		}catch(IOException ex){
			ex.printStackTrace();
		}
	}
	
	public void writeConfigToFile(String FileName){
		try {
			File pathFile = new File(FileName);
			DataOutputStream dos = new DataOutputStream(new FileOutputStream(pathFile, true));
			for (int i = 0; i < ising.Lp*ising.Lp; i ++){
				dos.writeInt(i);
				dos.writeChar('\t');
				dos.writeDouble(ising.phi[i]);
				dos.writeChar('\n');
			}
			dos.close();
		} catch (IOException ex){
			ex.printStackTrace();
		}
	}

	public void initFile(String file, boolean SvH){
		FileUtil.deleteFile(file);
		if(SvH){
			FileUtil.printlnToFile(file, " # SF vs H data ");			
			FileUtil.printlnToFile(file, " # Temperature = ", ising.T);
			FileUtil.printlnToFile(file, " # Data = H, S(k*), Free Energy, time");
		}else{
			FileUtil.printlnToFile(file, " # SF vs T data ");				
			//FileUtil.printlnToFile(file, " # External field = ", ising.H);
			FileUtil.printlnToFile(file, " # Data = H, S(k*), Free Energy, time");
		}
		FileUtil.printlnToFile(file, " # Density = ", ising.DENSITY);		
	}

}

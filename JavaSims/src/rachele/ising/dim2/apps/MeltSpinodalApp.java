package rachele.ising.dim2.apps;

import static java.lang.Math.PI;
import static java.lang.Math.floor;
import static java.lang.Math.sqrt;
import static scikit.numerics.Math2.j1;
import static scikit.util.Utilities.frame;
import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import rachele.ising.dim2.IsingField2D;
import rachele.ising.dim2.StructureFactor;
import rachele.util.FileUtil;
import scikit.graphics.dim2.Grid;
import scikit.jobs.Control;
import scikit.jobs.Job;
import scikit.jobs.params.ChoiceValue;
import scikit.jobs.params.DoubleValue;
import scikit.jobs.Simulation;
import static java.lang.Math.*;
//import rachele.util.FileUtil;

public class MeltSpinodalApp extends Simulation{
    Grid grid = new Grid("Phi(x)");
	StructureFactor sf;
    IsingField2D ising;
	boolean initFile = false;
    public int lastClear;

    
	public static void main(String[] args) {
		new Control(new MeltSpinodalApp(), "Ising Field");
	}
	
	public MeltSpinodalApp() {
		
		frame(grid);
		params.addm("Zoom", new ChoiceValue("Yes", "No"));
		params.addm("Interaction", new ChoiceValue( "Circle","Square"));
		params.addm("Noise", new ChoiceValue("Off","On"));
		params.addm("Dynamics?", new ChoiceValue("Langevin Conserve M", "Conjugate Gradient Min", 
				"Steepest Decent",  "Langevin No M Convervation"));
		params.add("Init Conditions", new ChoiceValue("Random Gaussian", 
				"Artificial Stripe 3", "Read From File", "Constant" ));
		params.addm("Approx", new ChoiceValue("Exact Stable",
				"Avoid Boundaries", "Exact SemiStable", "Exact", "Linear",  "Phi4"));
		params.addm("Plot FEvT", new ChoiceValue("Off", "On"));
		params.addm("Horizontal Slice", new DoubleValue(0.5, 0, 0.9999).withSlider());
		params.addm("Vertical Slice", new DoubleValue(0.5, 0, 0.9999).withSlider());
		params.addm("kR", new DoubleValue(5.135622302, 0.0, 6.0).withSlider());
		params.addm("T", 0.04);
		params.addm("H", 0.0);
		params.addm("dT", 0.001);
		params.addm("tolerance", 0.0001);
		params.addm("dt", 1.0);
		params.addm("J", -1.0);
		params.addm("R", 1000000.0);
		params.add("L/R", 6.0);
		params.add("R/dx", 16.0);
		params.add("kR bin-width", 0.1);
		params.add("Random seed", 0);
		params.addm("Magnetization", 0.75);
		params.add("time");
		flags.add("Clear");
	}
	
	public void animate() {
		params.set("time", ising.time());
		
		if (params.sget("Zoom").equals("Yes")) 
			grid.setAutoScale();
		else 
			grid.setScale(-1, 1);
		grid.registerData(ising.Lp, ising.Lp, ising.phi);
		
		if (flags.contains("Clear")){// || lastClear > 1000) {
			ising.getFreeEnergyAcc().clear();
			sf.getPeakH().clear();
			sf.getPeakV().clear();
			sf.getPeakC().clear();
			sf.getPeakHslope().clear();
			sf.getPeakVslope().clear();
			sf.getAccumulatorVA().clear();
			sf.getAccumulatorHA().clear();
			ising.aveCount = 0;
			lastClear = 0;
		}
		
		flags.clear();
	}
	
	public void clear() {
		initFile = false;
	}
	
	public void run() {
		ising = new IsingField2D(params);
		double binWidth = params.fget("kR bin-width");
		binWidth = IsingField2D.KR_SP / floor(IsingField2D.KR_SP/binWidth);
        sf = new StructureFactor(ising.Lp, ising.L, ising.R, binWidth, ising.dt);
		sf.setBounds(0.1, 14);
		double dT = .001;
		//double kR =5.135622302;
		String meltFile = "../../../research/javaData/sfData/Smelt";
		for (double den = -0.8; den < .95; den = den + .1){
			//double h=.2;
			params.set("Magnetization", den);
			ising.DENSITY = den;
			ising.randomizeField(0);
			double spinTemp=findOrderSpinodal(den);
			double temp=setTemp(spinTemp,dT);
			//Equilibrate			
			double lastSfValue = 0;
			for(int i = 0; i < 5000; i ++){
				ising.simulate();
				Job.animate();
//				if(ising.t%100 == 0){
//					sf.accumulateMelt(ising.phi);
//					writeDataToFile();
//				}
			}	
			boolean circleOn;
			if(params.sget("Interaction")=="Circle")
							circleOn=true;
			else
				circleOn=false;
			boolean vertStripes = true;
			int maxi=1;
			if(params.sget("Interaction")=="Square"){
				vertStripes = sf.findStripeDirection(ising.phi);
				if(vertStripes==true)
					System.out.println("Square interaction: vertical Stripes");
				else
					System.out.println("Square interaction: horizontal Stripes");
			}else{
				maxi = sf.clumpsOrStripes(ising.phi);
			}
			sf.accumulateMelt(circleOn,ising.phi,maxi);
			double sfValue = findSfValue(vertStripes);
			System.out.println(sfValue);
			while(sfValue > 100){
				for(int i = 0; i < 5000; i ++){
					ising.simulate();
					Job.animate();
					if(ising.time()%100 == 0){
						sf.accumulateMelt(circleOn,ising.phi,maxi);
						writeDataToFile();
					}
				}
				lastSfValue = sfValue;
				sf.accumulateMelt(circleOn,ising.phi,maxi);
				sfValue = findSfValue(vertStripes);
				while(abs(sfValue-lastSfValue) > (sfValue+lastSfValue)*.001){
					//while(abs(sfValue-lastSfValue) > (sfValue+lastSfValue)*(.1)){
					lastSfValue = sfValue;
					for(int i = 0; i < 1000; i ++){
						ising.simulate();
						params.set("time", ising.time());
						Job.animate();
					}
					sf.accumulateMelt(circleOn,ising.phi,maxi);
					writeDataToFile();
					sfValue = findSfValue(vertStripes);
					double diff = abs(sfValue-lastSfValue);
					System.out.println("diff = " + diff);
					Job.animate();
				}
				String equilFile="../../../research/javaData/sfData/equil";
				double kR = findkRvalue(maxi);
				FileUtil.printlnToFile(equilFile, ising.T, sf.peakValueC(), ising.freeEnergy, kR);
				temp += dT;
				ising.T = temp;
				System.out.println("temp = " + temp);
			}

			double meltTemp = temp -dT -dT/2;
			double error = dT/2;
			System.out.println("DONE:  Spinodal Temp  for density " +ising.DENSITY + " = "+meltTemp + " " + sfValue + " " + ising.t);
			double kR = findkRvalue(maxi);
			FileUtil.printlnToFile(meltFile, den, meltTemp, error, kR);
			}
					
        
 	}
	
	private double findkRvalue(int maxi){
		double kRvalue=0;
		for (int y = -ising.Lp/2; y < ising.Lp/2; y++) {
			for (int x = -ising.Lp/2; x < ising.Lp/2; x++) {
				double kR = (2*PI*sqrt(x*x+y*y)/ising.L)*ising.R;
				int i = ising.Lp*((y+ising.Lp)%ising.Lp) + (x+ising.Lp)%ising.Lp;
				if (i==maxi)
					kRvalue=kR;
				
			}
		}
		return kRvalue;
	}

	private double findSfValue(boolean veticalStripes){
		double sfvalue=0;
		if(params.sget("Interaction")=="Circle"){
			sfvalue = sf.peakValueC();	
		}else if(params.sget("Interaction")=="Square"){
			if(veticalStripes==true)
				sfvalue = sf.peakValueH();
			else
				sfvalue = sf.peakValueV();
		}		
		return sfvalue;
	}
	
	private double setTemp(double Ts, double dT) {
//		double sub = .003;
//		double temp = Ts-sub;
//		while(temp < 0){
//			sub /= 2;
//			temp = Ts-sub;
//		}
		System.out.println(Ts);
		double temp = (5.0/6.0)*Ts;
		temp /= dT;
		temp = (double)(int)temp;
		temp *= dT;
		ising.T = temp;
		params.set("T", temp);

		return temp;
	}

	private double findOrderSpinodal(double den) {
		double Ts=0;
		if(params.sget("Interaction")=="Circle"){
			double kR = 5.135622302;
			Ts = -circlePotential(kR)*(1-pow(den,2));
		}else if(params.sget("Interaction")=="Square"){
			double kR = 4.4934092;
			Ts = -squarePotential(kR)*(1-pow(den,2));	
		}
			return Ts;
	}
	private double circlePotential(double kR){
		return 2*j1(kR)/kR;
	}

	private double squarePotential(double kR){
		return sin(kR)/kR;
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
			params.set("H", readData);
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
	
	public void writeInputParams(String FileName){
		try {
			File inputFile = new File(FileName);
			DataOutputStream dos = new DataOutputStream(new FileOutputStream(inputFile, true));
			
			dos.writeDouble(params.fget("J"));
			dos.writeChar('\t');
			dos.writeDouble(params.fget("H"));
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

	public void initFile(String file, boolean SvH){
		FileUtil.deleteFile(file);
		if(SvH){
			FileUtil.printlnToFile(file, " # SF vs H data ");			
			FileUtil.printlnToFile(file, " # Temperature = ", ising.T);
			FileUtil.printlnToFile(file, " # Data = H, S(k*), Free Energy, time");
		}else{
			FileUtil.printlnToFile(file, " # SF vs T data ");				
			FileUtil.printlnToFile(file, " # External field = ", ising.H);
			FileUtil.printlnToFile(file, " # Data = H, S(k*), Free Energy, time");
		}
		FileUtil.printlnToFile(file, " # Density = ", ising.DENSITY);		
	}
	
	public void writeDataToFile(){
		boolean SvH = false;
		if (params.sget("Interaction")=="Square"){
				String dataFileV = "../../../research/javaData/sfData/dataFileV";
				String dataFileH = "../../../research/javaData/sfData/dataFileH";
				if (initFile == false){
					initFile(dataFileV, SvH);
					initFile(dataFileH, SvH);
					initFile = true;
				}
				if (SvH){
					FileUtil.printlnToFile(dataFileH, ising.H, sf.peakValueH(), ising.freeEnergy, ising.time());
					FileUtil.printlnToFile(dataFileV, ising.H, sf.peakValueV(), ising.freeEnergy, ising.time());					
				}else{
					FileUtil.printlnToFile(dataFileH, ising.T, sf.peakValueH(), ising.freeEnergy, ising.time());
					FileUtil.printlnToFile(dataFileV, ising.T, sf.peakValueV(), ising.freeEnergy, ising.time());
				}			
				System.out.println("Data written to file for time = " + ising.time());
		}else if(params.sget("Interaction")== "Circle"){
			String dataStripe = "../../../research/javaData/sfData/dataStripe";
			String dataClump = "../../../research/javaData/sfData/dataClump";
			if (initFile == false){
				initFile(dataStripe, SvH);
				initFile(dataClump, SvH);
				initFile = true;
			}
			if(SvH){
				FileUtil.printlnToFile(dataClump, ising.H, sf.peakValueC(), ising.freeEnergy, ising.time());
				FileUtil.printlnToFile(dataStripe, ising.H, sf.peakValueS(), ising.freeEnergy, ising.time());					
			}else{
				FileUtil.printlnToFile(dataClump, ising.T, sf.peakValueC(), ising.freeEnergy, ising.time());
				FileUtil.printlnToFile(dataStripe, ising.T, sf.peakValueS(), ising.freeEnergy, ising.time());
			}
			//System.out.println("Data written to file for time = " + ising.time());
		}
	}

}

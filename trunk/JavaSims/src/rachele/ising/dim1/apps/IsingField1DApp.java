package rachele.ising.dim1.apps;


import static java.lang.Math.PI;
import static java.lang.Math.floor;
import static scikit.util.Utilities.format;
//import static scikit.util.Utilities.frameTogether;
import scikit.dataset.PointSet;
import scikit.jobs.Control;
import scikit.jobs.Job;
import scikit.jobs.Simulation;
import scikit.jobs.params.ChoiceValue;
import scikit.jobs.params.DoubleValue;
import scikit.graphics.dim2.Plot;
import rachele.ising.dim1.FieldIsing1D;
import rachele.ising.dim1.StructureFactor1D;
import rachele.util.*;
import java.awt.Color;
import java.io.*;


public class IsingField1DApp extends Simulation{

	Plot fieldPlot = new Plot("Coarse Grained Field", "x", "phi(x)");
    Plot SFPlot = new Plot("Structure factor");
    Plot freeEngDenPlot = new Plot("Free Energy Density");
    Plot freeEngPlot = new Plot("Free Energy");
    Plot checkFT = new Plot("check FT");
    FieldIsing1D ising;
    StructureFactor1D sf;
    public int timeCount;
    
	public static void main(String[] args) {
		new Control(new IsingField1DApp(), "Ising Field 1D");
	}
	
	public void load(Control c) {
		//c.frameTogether("Displays", fieldPlot, freeEngPlot, freeEngDenPlot, SFPlot);
		c.frame(fieldPlot);
		c.frame(checkFT);
		
		////	Default parameters for nucleation
//		//params.addm("Model", new ChoiceValue("A", "B"));
//		params.addm("Noise", new ChoiceValue("On", "Off"));
//		params.addm("T", 0.85);
//		params.addm("J", -1.0);
//		params.addm("H", 0.04);
//		params.addm("R", 100000);
//		params.add("L/R", 32.0);
//		params.add("R/dx", 4.0);
//		params.add("kR bin-width", 0.1);
//		params.add("Random seed", 0);
//		params.add("Density", -.4);
//		params.add("dt", 0.1);
//		params.add("Time Allocation");
//		params.add("max Write Time", 30.0);
//		params.add("Time Count");
//		params.add("Time");
//		params.add("DENSITY");
//		params.add("Lp");
//		params.add("F");
		
//  Default params for clump model
		//params.addm("Model", new ChoiceValue("A", "B"));
		params.addm("Noise", new ChoiceValue("On", "Off"));
		params.addm("Interaction", new ChoiceValue("Step Function", "Effective Circle"));
		params.addm("T", new DoubleValue(0.04, 0, 0.2).withSlider());
		//params.addm("T", 0.04);
		params.addm("J", +1.0);
		params.addm("H", 0.8);
		params.addm("R", 2000000);
		params.add("L/R", 3.0);
		params.add("R/dx", 42.9);
		params.add("kR bin-width", 0.1);
		params.add("Density", -.4);
		params.add("Time Allocation");
		params.add("max Write Time", 30.0);
		params.add("Ry/Rx for Circle", 0.86747);
		params.add("Time");
		params.add("DENSITY");
		params.add("Lp");
		//params.add("Time Count");
		//params.add("F");		
		//flags.add("Write");
		flags.add("Write Config");
	}
	
	public void animate() {
		params.set("Time", format(ising.t));
		//params.set("F", format(ising.freeEnergyDensity));
		ising.readParams(params);
		//fieldPlot.setAutoScale(true);
		//fieldPlot.set
		
		checkFT.registerLines("checking", new PointSet(0,1,findCircleIntFunction(params.fget("Ry/Rx for Circle"))), Color.BLUE);
		SFPlot.registerLines("Structure factor", sf.getAccumulator(), Color.BLACK);
		fieldPlot.registerLines("Field", new PointSet(0, ising.dx, ising.phi), Color.BLACK);
		freeEngDenPlot.registerLines("F.E. density", ising.getFreeEngAcc(), Color.BLACK);
		freeEngPlot.registerLines("Free energy", new PointSet(0, ising.dx, ising.F), Color.BLACK);
		
		if (flags.contains("Clear S.F.")) {
			sf.getAccumulator().clear();
			System.out.println("clicked");
		}
		flags.clear();
	}
	
	public void clear() {
		fieldPlot.clear();
		SFPlot.clear();
		freeEngDenPlot.clear();
		freeEngPlot.clear();
	}
	
	public void run(){
		ising = new FieldIsing1D(params);
		boolean circleInt = false;
		if(params.sget("Interaction")=="Effective Circle")
			circleInt = true;
		double [] circleIntFT = new double [ising.Lp];
		if (circleInt) circleIntFT = findCircleIntFunction(params.fget("Ry/Rx for Circle"));
		//String pathFileName = "racheleDataFiles/inputPath";
		//String inputFileName = "racheleDataFiles/inputParams";
		String inputFileName = "../../../research/javaData/configs1d/inputConfig";
		String configFileName = "../../../research/javaData/configs1d/config";
	    int maxWriteCount = (int)(params.fget("max Write Time")/ising.dt);
	    params.set("Time Allocation", maxWriteCount);
	    double maxWriteTime = maxWriteCount*ising.dt;
	    params.set("max Write Time", maxWriteTime);
		double KR_SP = FieldIsing1D.KR_SP;
		double binWidth = KR_SP / floor(KR_SP/params.fget("kR bin-width"));

		sf = new StructureFactor1D(ising.Lp, ising.L, ising.R, binWidth);
		Job.animate();
		
		//sf.getAccumulator().clear();
		timeCount = maxWriteCount +1;
		
		while (true) {
			if (flags.contains("Write")) {
				timeCount = 0;
				FileUtil.deleteFile(configFileName);
				FileUtil.deleteFile(inputFileName);
				writeInputParams(inputFileName);
				while (timeCount <= maxWriteCount){
					if(circleInt)
						ising.simulateCircle(circleIntFT);
					else
						ising.simulate();			
					writeConfigToFileWithTime(configFileName);				
					Job.animate();					
					params.set("Time Count", timeCount);
					timeCount += 1;
				}
				flags.clear();
			}else if(flags.contains("Write Config")){
				FileUtil.deleteFile(configFileName);
				FileUtil.writeConfigToFile(configFileName, ising.Lp, ising.phi);
			}
			
			ising.simulate();
			sf.accumulate(ising.phi);
			Job.animate();
			//SFPlot.setDataSet(0, sf.getAccumulator());
			//fieldPlot.setDataSet(0, new PointSet(0, ising.dx, ising.phi));
			
		}
	}
	
	/**
	 * Find a value for the FT for each allowed kR: 
	 * kR = 2 PI i R / L
	 */
	public double [] findCircleIntFunction(double ampFactor){
		double [] circleFT = new double [ising.Lp];
		int N = 3;
		double minx = -1.0; double maxx = 1.0;
		double range = maxx - minx;
		double delta = range/(double)N;
		for (int i = 0; i < ising.Lp; i++){
			double kRvalue = (2*PI*i/ising.L) * ising.R;
			double coskR = Math.cos(kRvalue);
			// have to intergrate: int _-1 ^1 cos(kR*x)*sqrt(1-x*x)
			double sum = 0;
			for(int j = 0; j < N; j++){
				double xValue = (double)j*delta-1.0;
				double factor = Math.sqrt(1-xValue*xValue);
				sum += coskR*factor;
				//System.out.println(i + " " +coskR+ " " + factor + " " + sum);
			}
			sum*= delta;
			circleFT[i] = sum*ampFactor;
//			System.out.println(i + " " +circleFT[i]+ " " + );
		}
		return circleFT;
	}
	
	public void writeInputParams(String FileName){
		try {
			File inputFile = new File(FileName);
			DataOutputStream dos = new DataOutputStream(new FileOutputStream(inputFile, true));
			//System.out.println(params.fget("T"));
			dos.writeDouble(params.fget("T"));
			dos.writeChar('\t');
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
			dos.writeDouble(params.fget("dt"));
			dos.writeChar('\t');
			dos.writeDouble(params.fget("max Write Time"));
			dos.writeChar('\n');
			dos.close();
		}catch(IOException ex){
			ex.printStackTrace();
		}
	}

	public void writeConfigToFileWithTime(String FileName){
		try {
			File pathFile = new File(FileName);
			DataOutputStream dos = new DataOutputStream(new FileOutputStream(pathFile, true));
			for (int i = 0; i < ising.Lp; i ++){
				dos.writeInt(timeCount);
				dos.writeChar('\t');
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
	
}

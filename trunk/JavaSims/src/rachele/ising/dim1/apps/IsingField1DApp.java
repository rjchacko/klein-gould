package rachele.ising.dim1.apps;

import static java.lang.Math.floor;
import static scikit.util.Utilities.format;
import scikit.dataset.PointSet;
import scikit.jobs.Control;
import scikit.jobs.Job;
import scikit.jobs.Simulation;
import scikit.jobs.params.ChoiceValue;
import scikit.jobs.params.DirectoryValue;
import scikit.jobs.params.DoubleValue;
import scikit.graphics.dim2.Plot;
import rachele.ising.dim1.FieldIsing1D;
import rachele.ising.dim1.StructureFactor1D;
import rachele.util.*;
import java.awt.Color;
import java.io.*;

public class IsingField1DApp extends Simulation{

	Plot fieldPlot = new Plot("Coarse Grained Field", "x", "phi(x)");
	Plot ftPlot = new Plot("FT", "k", "ft_phi(x)");
    Plot SFPlot = new Plot("Structure factor");
    Plot test = new Plot("test phi", "x", "phi(x)");
    Plot testPlot = new Plot("test phi ft","k", "ft_phi(x)");
    FieldIsing1D ising;
    StructureFactor1D sf;
    public int timeCount;
    
	public static void main(String[] args) {
		new Control(new IsingField1DApp(), "Ising Field 1D");
	}
	
	public void load(Control c) {
		//c.frameTogether("Displays", fieldPlot, freeEngPlot, freeEngDenPlot, SFPlot);
		c.frameTogether("data", fieldPlot, ftPlot, test, testPlot);
		//	Default parameters for nucleation
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

		params.add("Config Directory",new DirectoryValue("/home/erdomi/data/lraim/configs1dAutoName"));
		params.addm("Noise", new DoubleValue(1.0, 0, 1.0).withSlider());
		params.addm("Dynamics", new ChoiceValue("Langevin","Conserved F Space", "Conserved","Conserved w mob", "Conserved semi imp",  "Glauber"));
		params.addm("Random Seed", 0);
		params.addm("T", new DoubleValue(0.08, 0, 0.2).withSlider());
		params.addm("J", +1.0);
		params.addm("H", 0.65);
		params.addm("R", 4600000);
		params.add("L/R", 2.7826087);
		params.add("R/dx", 46.0);
		params.add("kR bin-width", 0.1);
		params.add("Density", -.4);
		params.add("Time Allocation");
		params.add("max Write Time", 30.0);
		params.addm("dt",0.1);
		params.add("Time");
		params.add("DENSITY");
		params.add("Lp");
		params.add("Free Energy");
		
		flags.add("Write Config");
	}
	
	public void animate() {
		params.set("Time", format(ising.t));
		params.set("Free Energy", format(ising.freeEnergyDensity));
		ising.readParams(params);
		fieldPlot.setAutoScale(true);
		test.setAutoScale(true);
		testPlot.setAutoScale(true);
		
		SFPlot.registerLines("Structure factor", sf.getAccumulator(), Color.BLACK);
		fieldPlot.registerLines("Field", new PointSet(0, ising.dx, ising.phi), Color.BLACK);
		double [] realFT = new double [ising.Lp];
		double [] imagFT = new double [ising.Lp];
		for (int i = 1; i < ising.Lp; i++){
			realFT[i] = ising.phi_k[2*i];
			imagFT[i] = ising.phi_k[2*i+1];
		}
		
		ftPlot.registerLines("real", new PointSet(0, ising.dx, realFT), Color.BLUE);	
		ftPlot.registerLines("imag", new PointSet(0, ising.dx, imagFT), Color.RED);	
		
		for (int i = 1; i < ising.Lp; i++){
			realFT[i] = ising.phi2_k[2*i];
			imagFT[i] = ising.phi2_k[2*i+1];
		}
		testPlot.registerLines("real2", new PointSet(0, ising.dx, realFT), Color.BLUE);	
		testPlot.registerLines("imag2", new PointSet(0, ising.dx, imagFT), Color.RED);			
		test.registerLines("tester", new PointSet(0, ising.dx, ising.phi2), Color.RED);
//		testPlot.registerLines("test ft", new PointSet(0, ising.dx, ising.phi2_k), Color.GREEN);
		
		if (flags.contains("Clear S.F.")) {
			sf.getAccumulator().clear();
			System.out.println("clicked");
		}
		flags.clear();
	}
	
	public void clear() {
		SFPlot.clear();
	}
	
	public void run(){
		ising = new FieldIsing1D(params);
		String writeDir = params.sget("Config Directory");
		StringBuffer sb = new StringBuffer();
		sb.append(writeDir); sb.append(File.separator); sb.append("L"); sb.append(ising.Lp);
		int range = (int)(ising.Lp/params.fget("L/R")); sb.append("R"); sb.append(range);
		sb.append("T"); sb.append(ising.T); sb.append("h"); sb.append(ising.H);
		String configFileName = sb.toString();
	    sb.append("params"); String paramsFile = sb.toString();
		int maxWriteCount = (int)(params.fget("max Write Time")/ising.dt);
	    params.set("Time Allocation", maxWriteCount);
	    double maxWriteTime = maxWriteCount*ising.dt;
	    params.set("max Write Time", maxWriteTime);
		double KR_SP = FieldIsing1D.KR_SP;
		double binWidth = KR_SP / floor(KR_SP/params.fget("kR bin-width"));

		sf = new StructureFactor1D(ising.Lp, ising.L, ising.R, binWidth);
		Job.animate();
		
		timeCount = maxWriteCount +1;
		boolean glauber;
		if(params.sget("Dynamics")=="Glauber")
			 glauber = true;
		else
			glauber = false;
		System.out.println("Galuber dynamics is " + glauber);
		

		
		while (true) {
			if (flags.contains("Write")) {
				timeCount = 0;
				FileUtil.deleteFile(configFileName);
				FileUtil.initFile(paramsFile, params);
				while (timeCount <= maxWriteCount){
					writeConfigToFileWithTime(configFileName);				
					Job.animate();					
					params.set("Time Count", timeCount);
					timeCount += 1;
				}
				flags.clear();
			}else if(flags.contains("Write Config")){
				FileUtil.initFile(paramsFile, params);
				FileUtil.deleteFile(configFileName);
				FileUtil.writeConfigToFile(configFileName, ising.Lp, ising.phi);
			}
			if(params.sget("Dynamics")=="Glauber") ising.simulateGlauber();
			else if(params.sget("Dynamics")=="Langevin")ising.simulate();
			else if (params.sget("Dynamics")=="Conserved"){ ising.simulateConserved();
//				ising.simulateConservedFspace();
//				ising.simulateConseveredSemiImp();
			}else if (params.sget("Dynamics")=="Conserved w mob") ising.simulateConservedWithMobility();
			else if (params.sget("Dynamics")=="Conserved F Space"){
				ising.simulateConserved();
				ising.simulateConservedFspace();
			}
			else if (params.sget("Dynamics")=="Conserved semi imp") ising.simulateConseveredSemiImp();			
			sf.accumulate(ising.phi);
			Job.animate();
			//SFPlot.setDataSet(0, sf.getAccumulator());
			//fieldPlot.setDataSet(0, new PointSet(0, ising.dx, ising.phi));
			
		}
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

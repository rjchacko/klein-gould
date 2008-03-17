package rachele.ising.dim1.apps;


import static java.lang.Math.floor;
import static scikit.util.Utilities.format;
import static scikit.util.Utilities.frameTogether;
import scikit.dataset.PointSet;
import scikit.jobs.Control;
import scikit.jobs.Job;
import scikit.jobs.Simulation;
import scikit.jobs.params.ChoiceValue;
import scikit.graphics.dim2.Plot;
import rachele.ising.dim1.FieldIsing1D;
import rachele.ising.dim1.StructureFactor1D;

import java.awt.Color;
import java.io.*;


public class IsingField1DApp extends Simulation{

	Plot fieldPlot = new Plot("Coarse Grained Field");
    Plot SFPlot = new Plot("Structure factor");
    Plot freeEngDenPlot = new Plot("Free Energy Density");
    Plot freeEngPlot = new Plot("Free Energy");
    FieldIsing1D ising;
    StructureFactor1D sf;
    public int timeCount;
    
	public static void main(String[] args) {
		new Control(new IsingField1DApp(), "Ising Field 1D");
	}
	
	public IsingField1DApp(){
		frameTogether("Displays", fieldPlot, freeEngPlot, freeEngDenPlot, SFPlot);
//	Defoult parameters for nucleation
		params.addm("Zoom", new ChoiceValue("A", "B"));
		params.addm("Noise", new ChoiceValue("On", "Off"));
		params.addm("T", 0.85);
		params.addm("J", -1.0);
		params.addm("H", 0.04);
		params.addm("R", 100000);
		params.add("L/R", 32.0);
		params.add("R/dx", 4.0);
		params.add("kR bin-width", 0.1);
		params.add("Random seed", 0);
		params.add("Density", -.4);
		params.add("dt", 0.1);
		params.add("Time Allocation");
		params.add("max Write Time", 30.0);
		params.add("Time Count");
		params.add("Time");
		params.add("DENSITY");
		params.add("Lp");
		params.add("F");
		
		flags.add("Write");

//Default Parameters for clumps
//		params.addm("Zoom", new ChoiceValue("B", "A"));
//		params.addm("T", 0.1);
//		params.addm("J", 1.0);
//		params.addm("dt", 0.1);
//		params.addm("R", 3000);
//		params.addm("H", 0.00);
//		params.add("L/R", 16.0);
//		params.add("R/dx", 16.0);
//		params.add("kR bin-width", 0.1);
//		params.add("Random seed", 0);
//		params.add("Density", 0.0);
//		params.add("Time");
//		params.add("DENSITY");
//		params.add("Lp");	
//		params.add("F");
		
	}
	
	public void animate() {
		params.set("Time", format(ising.t));
		params.set("F", format(ising.freeEnergyDensity));
		ising.readParams(params);
		
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
	
	public void deleteFile(String fileName){
		File file = new File(fileName);
		boolean success = file.delete();
		if (success)
			System.out.println("File deleted");
		else
			System.out.println("File delete failed");			
	}

	public void writeConfigToFile(String FileName){
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

	public void run(){
		ising = new FieldIsing1D(params);
		
		String pathFileName = "racheleDataFiles/inputPath";
		String inputFileName = "racheleDataFiles/inputParams";
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
				deleteFile(pathFileName);
				deleteFile(inputFileName);
				writeInputParams(inputFileName);
				while (timeCount <= maxWriteCount){
					ising.simulate();			
					writeConfigToFile(pathFileName);				
					Job.animate();					
					params.set("Time Count", timeCount);
					timeCount += 1;
				}
				flags.clear();
			}
			ising.simulate();
			sf.accumulate(ising.phi);
			Job.animate();
			//SFPlot.setDataSet(0, sf.getAccumulator());
			//fieldPlot.setDataSet(0, new PointSet(0, ising.dx, ising.phi));
			
		}
	}
	
}

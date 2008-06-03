package rachele.ising.dim1.apps;

import static scikit.util.Utilities.format;
import scikit.dataset.PointSet;
import scikit.jobs.Control;
import scikit.jobs.Job;
import scikit.jobs.Simulation;
import scikit.jobs.params.ChoiceValue;
import scikit.jobs.params.DoubleValue;
import scikit.graphics.dim2.Plot;
import rachele.ising.dim1.FieldIsing1D;
import rachele.util.*;
import static scikit.util.DoubleArray.*;
import java.awt.Color;
import java.io.*;

public class IsingField1DCircleApp extends Simulation{


	Plot fieldPlot = new Plot("Coarse Grained Field", "x", "phi(x)");
	FieldIsing1D ising;
	public int timeCount;
    
	public static void main(String[] args) {
		new Control(new IsingField1DCircleApp(), "Effective Circle Interaction Simulation");
	}
	
	public void load(Control c) {
		c.frame(fieldPlot);
			
		params.addm("Noise", new ChoiceValue("On", "Off"));
		params.addm("Random Seed", 0);
		params.addm("T", new DoubleValue(0.06, 0, 0.2).withSlider());
		params.addm("J", +1.0);
		params.addm("H", 0.81);
		params.addm("Target Mean Phi", 0.435);
		params.addm("R", 2490000.0);
		params.add("L/R", 2.409638554);
		params.add("R/dx", 80.0);
		params.add("kR bin-width", 0.1);
		params.add("Density", -.4);
		params.add("Time Allocation");
		params.add("max Write Time", 30.0);
		params.addm("Ry/Rx for Circle", 1.0);
		params.add("Time");
		params.add("DENSITY");
		params.add("Lp");
		
		flags.add("Write Config");
	}
	
	public void animate() {
		params.set("Time", format(ising.t));
		ising.readParams(params);
		fieldPlot.registerLines("Field", new PointSet(0, ising.dx, ising.phi), Color.BLACK);
	}
	
	public void clear() {
		fieldPlot.clear();
	}
	
	public void run(){
		ising = new FieldIsing1D(params);
		double ampFactor=1;
		ampFactor = params.fget("Ry/Rx for Circle"); 			
		String inputFileName = "../../../research/javaData/configs1d/inputConfigOpt";
		String configFileName = "../../../research/javaData/configs1d/configOpt";

		Job.animate();

		while (true) {
			if (flags.contains("Write")) {
				timeCount = 0;
				FileUtil.deleteFile(configFileName);
				FileUtil.deleteFile(inputFileName);
				writeInputParams(inputFileName);
			}else if(flags.contains("Write Config")){
				FileUtil.deleteFile(configFileName);
				FileUtil.writeConfigToFile(configFileName, ising.Lp, ising.phi);
			}			
			flags.clear();
			ampFactor = params.fget("Ry/Rx for Circle"); 
			ising.simulateCircle(ampFactor);
			double newH = 0;
			double deltaH = 0.00001;
			if (mean(ising.phi) < params.fget("Target Mean Phi"))
				newH = ising.H + deltaH;
			else
				newH = ising.H - deltaH;
			params.set("H", newH);
			Job.animate();			
		}
	}
		
	public void writeInputParams(String FileName){
		try {
			File inputFile = new File(FileName);
			DataOutputStream dos = new DataOutputStream(new FileOutputStream(inputFile, true));
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

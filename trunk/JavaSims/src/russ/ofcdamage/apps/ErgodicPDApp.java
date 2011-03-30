package russ.ofcdamage.apps;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.text.DecimalFormat;

import scikit.jobs.Control;
import scikit.jobs.Job;
import scikit.jobs.Simulation;
import scikit.jobs.params.DirectoryValue;
import chris.ofcdamage.basicOFC;
import chris.util.PrintUtil;

public class ErgodicPDApp extends Simulation {
	
	private static double chiTOL = 800;	

	private basicOFC model;
	private int pt, ptmax;
	private double aMin, time[], Metric[], dAA;
	private String dir, bname, fitfile;
	//private FitUtil fitter;
	private DecimalFormat fmt = new DecimalFormat("0.000");

	
	public static void main(String[] args) {
		new Control(new ErgodicPDApp(), "OFC Parameters");
	}
	
	
	public void load(Control c) {
		params.add("Data Directory",new DirectoryValue("/Users/cserino/Research/Data2/ErgPhaseDiagram/"));
		params.add("File Name","File_Name");
		params.add("Random Seed",0);
		params.add("Min Radius", (int) 1);
		params.add("Max Radius", (int) 5);
		params.add("Radius Step", (int) 1);
		params.add("Interaction Shape");
		params.add("Interaction Radius (R)");
		params.add("Lattice Size");
		params.add("Boundary Condtions");
		params.add("Failure Stress (\u03C3_f)");
		params.add("Residual Stress (\u03C3_r)");
		params.add("\u03C3_r width",0.21);
		params.add("\u03C3_r Step Size",0.01);
		params.add("Dissipation (\u03B1)");
		params.add("Equil Time");
		params.add("Trend Time");
		params.add("Number of Plate Updates");
		params.add("N_dead");

		params.set("Interaction Radius (R)",(int)(1));
		params.set("Interaction Shape","Circle");
		params.set("Boundary Condtions","Periodic");
		params.set("Lattice Size",1<<8);
		params.set("Equil Time", (int)(2e5));	// 2e5
		params.set("Trend Time", (int)(2e5));	// 2e5
		params.set("Failure Stress (\u03C3_f)",2.0);
		params.set("Residual Stress (\u03C3_r)",1.25);
		params.set("Dissipation (\u03B1)",0.01);

	}

	public void run() {

		int rmin = params.iget("Min Radius");
		int rmax = params.iget("Max Radius");
		int dRR  = params.iget("Radius Step");
		aMin     = params.fget("\u03C3_r width");
		dAA      = params.fget("\u03C3_r Step Size");
		dir      = params.sget("Data Directory");
		bname    = params.sget("File Name");
		fitfile  = dir + File.separator + bname + "_FitData.txt";
		
		PrintUtil.printlnToFile(fitfile,"Range","Noise","Chi^2");
			
			
		for (int rr = rmin ; rr <= rmax; rr += dRR ){
			params.set("Interaction Radius (R)",rr); 
			params.set("\u03C3_r width",aMin);
			getXover();
		}

		
		return;
	}
	
	
	private void getXover(){
		
		double aNOW = aMin;
		boolean ergodic = true;

		while(ergodic){
			ergodic = simSystem();
			params.set("Random Seed",(params.iget("Random Seed")+1));
			// Fix This (check)
			//aNOW += dAA;
			aNOW -= dAA;
			params.set("\u03C3_r width",aNOW);
		}
		
		aMin = aNOW + 3*dAA;
		
		return;
	}
	
	private boolean simSystem(){

		// Setup model
		params.set("N_dead","Initializing");
		model = new basicOFC(params);
		
		// Equilibrate
		equil();
		// Simulate w/o Damage for Data
		ideal();
		
		params.set("N_dead","Fitting");
		Job.animate();
		
		saveMetricData();
		// Fit and check for erdodicity
//		Metric = CopyUtil.invertAndScaleArray(Metric,Metric.length,1);
//		
//		fitter = new FitUtil(Metric.length);
//		
//		double[] foo = fitter.fit(time,Metric,1,false);
//		saveFitData(foo);
		
		// Fix This (check)
		return (0 < chiTOL);
		//return true;
	}
	
	
	private void saveMetricData(){
		
		int N = Metric.length;
		
		String fn = dir + File.separator + bname + "_" + params.iget("Interaction Radius (R)") + "_" + 
					fmt.format(params.fget("\u03C3_r width")) + ".txt";
		
		try{
			File file = new File(fn);
			PrintWriter pw = new PrintWriter(new FileWriter(file, true), true);
			pw.print("Time");
			pw.print("\t");
			pw.print("Metric");
			pw.println();
			for (int ii = 0 ; ii < N ; ii++){				
				pw.print(time[ii]);
				pw.print("\t");
				pw.print(Metric[ii]);
				pw.println();
			}
			pw.close();
		}
		catch (IOException ex){
			ex.printStackTrace();
		}
		
		return;
	}
	
	@SuppressWarnings("unused")
	private void saveFitData(double[] fitparams){

		PrintUtil.printlnToFile(fitfile,params.iget("Interaction Radius (R)"),params.fget("\u03C3_r width"),fitparams[4]);
		return;
	}
	
	private void equil(){
	
		// Equilibrate
		params.set("N_dead", "Equilibrating");
		pt      = 0;
		ptmax   = params.iget("Equil Time");
		
		model.runClocks(false);
		while(pt < ptmax){
			model.equilibrate();
			pt++;
			if(pt%200==0) Job.animate();
		}
		
		return;
	}
	
	private void ideal(){

		// Simulate w/o Damage for Data
		params.set("N_dead","Simulating w/o damage");
		model.runClocks(true);		
		int pt0 = 0;
		ptmax = params.iget("Trend Time"); 
		pt    = params.iget("Trend Time");
		
		time   = new double[ptmax];
		Metric = new double[ptmax]; 
		
		while(pt0 < ptmax){
			model.avalanche();
			time[pt0]     = model.getSimTime();
			Metric[pt0++] = model.getMetric();
			pt++;
			if(pt0%200==0) Job.animate();
		}
		
		return;
	}
	
	public void animate() {
		params.set("Number of Plate Updates", pt - ptmax);
		return;
	}

	public void clear() {
		return;
	}
	


}

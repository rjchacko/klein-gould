package chris.tests;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.text.DecimalFormat;

import scikit.jobs.Control;
import scikit.jobs.Job;
import scikit.jobs.Simulation;
import scikit.jobs.params.ChoiceValue;
import scikit.jobs.params.DirectoryValue;
import chris.MD.TwoD.LennardJones;
import chris.util.PrintUtil;

public class SFtestApp extends Simulation{
	
	public LennardJones model;
	static DecimalFormat tf = new DecimalFormat("######.00");
	static DecimalFormat Ef = new DecimalFormat("0.###E0");
	private double now, T, L, then, tau;
	private String fout;
	
	
	public static void main(String[] args) {
		new Control(new SFtestApp(), "Lennard-Jones System");
	}

	public void animate() {
		
		double[] eng = new double[2];
		if(T != params.fget("T")){
			T = params.fget("T");
			model.changeT(T);
		}
		tau = params.fget("d\u03C4");
		params.set("t", tf.format(now));
		model.getEsys(eng);
		params.set("T_m", Ef.format(eng[0]/model.N));
		params.set("E",Ef.format(eng[0]+eng[1]));	}

	public void clear() {
		
	}

	public void load(Control c) {
		
		params.add("Data Directory",new DirectoryValue("/Users/cserino/Desktop/"));
		params.add("Data File", "default");
		params.add("seed",0); 
		params.add("L",13.);
		params.add("Boundary Conditions", new ChoiceValue("Periodic","Closed"));
		params.add("Initial Conditions", new ChoiceValue("Debug","Read In", "Melt", "Copy", "Viscous"));
		params.add("ODE Solver", new ChoiceValue("Velocity Verlet","First Order Euler"));
		params.add("N",0);
		params.add("M",1);
		params.add("R",0.25); 
		params.add("dt",5e-3);
		params.addm("d\u03C4",50);
		params.addm("T",1e-6);
		params.add("t");
		params.add("E");
		params.add("T_m");
		
	}

//	public void run() {
//		model = new LennardJones(params);
//		L     = params.iget("L");
//		fout  = params.sget("Data Directory") + File.separator + params.sget("Data File") + ".txt";
//		int nmx = (int)(1e3);
//		double sf[][][]   = new double[1][nmx][nmx];
//
//		Job.animate();
//		
//		for(int nx = 0; nx < nmx ; nx++){	// using k^2 < 200
//			for(int ny = 0; ny < nmx ; ny++){	// using k^2 < 200
//				sf[0][nx][ny] = model.structureFactor(nx, ny, L, L);
//				params.set("t",ny);
//				Job.animate();
//			}
//			params.set("T_m",100.*nx/nmx);
//			Job.animate();
//		};
//
//		printSF(sf,nmx,1);
//		return;
//	}
	
	public void run() {
		model = new LennardJones(params);
		tau   = params.fget("d\u03C4");
		T     = params.fget("T");
		then  = 0;
		L     = params.fget("L");
		fout  = params.sget("Data Directory") + File.separator + params.sget("Data File") + ".txt";
		PrintUtil.printlnToFile(params.sget("Data Directory") + File.separator + "Params_" + params.sget("Data File") + ".log", params.toString());
		Job.animate();
		
		int tmx = 1000;
		int count = 0;
		int nmx = (int)(1e2);
		double sf[][][]   = new double[tmx][nmx][nmx];

		while(count < tmx){
			now = model.stepWT();
			Job.animate();
			if(now - then >= tau){
				then = now;
				for(int nx = 0; nx < nmx ; nx++){	// using k^2 < 200
					for(int ny = 0; ny < nmx ; ny++){	// using k^2 < 200
						sf[count][nx][ny] = model.structureFactor(nx, ny, L, L);
						params.set("t",ny);
						Job.animate();
					}
					params.set("T_m",100.*nx/nmx);
					Job.animate();
				}
				count++;
				Job.animate();
			}
		}

		printSF(sf,nmx,tmx);
		return;
	}
	
	private void printSF(double[][][]sf, int nmx, int tmx){
		try{
			File file = new File(fout);
			PrintWriter pw = new PrintWriter(new FileWriter(file, true), true);
			pw.println("t \t nx \t ny \t S(k)");
			for(int t = 0; t < tmx ; t++){
				for(int nx = 0; nx < nmx ; nx++){
					for(int ny = 0; ny < nmx ; ny++)
						pw.println(t+"\t"+nx+"\t"+ny+"\t"+sf[t][nx][ny]);
				}
			}
			pw.close();
		}
		catch (IOException ex){
			ex.printStackTrace();
		}
	return;
	}
}

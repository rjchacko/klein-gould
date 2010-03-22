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
	private double now, then, tau, T, L;
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
		params.add("Initial Conditions", new ChoiceValue("Read In", "Melt", "Copy", "Viscous"));
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
//
//		double sf[]   = new double[(int)(1e4)];
//
//		Job.animate();
//		
//		for(int nx = 0; nx < 1e4 ; nx++){	// using k^2 < 200
//			sf[nx] = model.structureFactor(nx, 0, L, L);
//			params.set("T_m",nx/1e2);
//			Job.animate();
//		};
//
//
//		try{
//			File file = new File(fout);
//			PrintWriter pw = new PrintWriter(new FileWriter(file, true), true);
//			pw.println("k \t S[k]");
//
//			for(int jj = 0 ; jj < 1e4 ; jj++){
//				pw.print(jj*Math.PI+"\t"); // L = 1;
//				pw.println(sf[jj]);
//			}
//			pw.close();
//		}
//		catch (IOException ex){
//			ex.printStackTrace();
//		}
//
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
		
		int count       = 0;
		int nkm         = 1000;
		int mxdp        = 100;
		double sf[][]   = new double[nkm][mxdp];
		int ny = 0;

		while(count < mxdp){
			now = model.stepWT();
			Job.animate();
			if(now - then >= tau){
				then = now;
				for(int nx = 0; nx < nkm ; nx++){
					sf[nx][count] = model.structureFactor(nx, ny, L, L);
				}
				count++;
				Job.animate();
			}
		}

		printSF(sf, mxdp);
		return;
	}
	
	private void printSF(double[][]sf, int lngth){
		try{
			File file = new File(fout);
			PrintWriter pw = new PrintWriter(new FileWriter(file, true), true);
			pw.println("nx \t S(k,t=0) \t S(k,t=1) \t S(k,t=2) \t . . . ");
			for(int nx = 0; nx < sf.length ; nx++){
				for(int cc = 0 ; cc < lngth ; cc++){
				pw.print(nx*Math.PI+"\t");
				pw.print(sf[nx][cc]+"\t");
				}
				pw.println();
			}
			pw.close();
		}
		catch (IOException ex){
			ex.printStackTrace();
		}
	return;
	}
}

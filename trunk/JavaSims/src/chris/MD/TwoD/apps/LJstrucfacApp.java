package chris.MD.TwoD.apps;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.text.DecimalFormat;

import chris.MD.TwoD.LennardJones;
import chris.util.PrintUtil;
import scikit.jobs.Control;
import scikit.jobs.Job;
import scikit.jobs.Simulation;
import scikit.jobs.params.ChoiceValue;
import scikit.jobs.params.DirectoryValue;

public class LJstrucfacApp extends Simulation{
	
	public LennardJones model;
	static DecimalFormat tf = new DecimalFormat("######.00");
	static DecimalFormat Ef = new DecimalFormat("0.###E0");
	private double now, then, tau, T, L;
	private String fout;
	
	
	public static void main(String[] args) {
		new Control(new LJstrucfacApp(), "Lennard-Jones System");
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
		params.set("E",Ef.format(eng[0]+eng[1]));
	}

	public void clear() {}

	public void load(Control c) {
		
		params.add("Data Directory",new DirectoryValue("/Users/cserino/Desktop/"));
		params.add("Data File", "default");
		params.add("seed",0); 
		params.add("L",15.);
		params.add("Boundary Conditions", new ChoiceValue("Periodic","Closed"));
		params.add("Initial Conditions", new ChoiceValue("Read In", "Melt", "Copy", "Viscous"));
		params.add("ODE Solver", new ChoiceValue("Velocity Verlet","First Order Euler"));
		params.add("N",0);
		params.add("M",1);
		params.add("R",0.25); 
		params.add("dt",5e-3);
		params.addm("d\u03C4",10);
		params.addm("T",1e-6);
		params.add("t");
		params.add("E");
		params.add("T_m");
		
	}

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
		int mxdp        = 100;
		int nkmx        = 1000;
		int nkmy        = nkmx;
		double sf[][][] = new double[nkmx][nkmy][mxdp];

		while(count < mxdp){
			now = model.stepWT();
			Job.animate();
			if(now - then >= tau){
				then = now;
				for(int nx = 0; nx < nkmx ; nx++){
					for(int ny = 0; ny < nkmy ; ny++){
						sf[nx][ny][count] = model.structureFactor(nx, ny, L, L);
					}
				}
				count++;
			}
		}

		printSF(sf, mxdp, nkmx, nkmy);
		return;
	}
	
	private void printSF(double[][][] sf, int lngth, int xm, int ym){
		try{
			File file = new File(fout);
			PrintWriter pw = new PrintWriter(new FileWriter(file, true), true);
			pw.println("t \t nx \t ny \t S(k)");
			for (int jj = 0 ; jj < lngth ; jj++){
				for(int nx = 0; nx < xm ; nx++){
					for(int ny = 0; ny < ym ; ny++){
						pw.print(jj+"\t");
						pw.print(nx+"\t");
						pw.print(ny+"\t");
						pw.println(sf[nx][ny][jj]);
					}
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

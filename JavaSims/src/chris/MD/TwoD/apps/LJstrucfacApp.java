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
	private double now, then, tau, T;
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
//		params.add("Ly",15.);
		params.add("Boundary Conditions", new ChoiceValue("Periodic","Closed"));
		params.add("Initial Conditions", new ChoiceValue("Read In", "Melt", "Copy", "Viscous","Debug"));
		params.add("ODE Solver", new ChoiceValue("Velocity Verlet","First Order Euler"));
		params.add("N",0);
		params.add("M",1);
		params.add("R",0.25); 
		params.add("dt",5e-3);
		params.addm("d\u03C4",5.);
		params.addm("T",1e-6);
		params.add("t");
		params.add("E");
		params.add("T_m");
		params.add("SF");
		
	}

	public void run() {
		model = new LennardJones(params);
		tau   = params.fget("d\u03C4");
		T     = params.fget("T");
		then  = 0;
		fout  = params.sget("Data Directory") + File.separator + params.sget("Data File") + ".txt";
		PrintUtil.printlnToFile(params.sget("Data Directory") + File.separator + "Params_" + params.sget("Data File") + ".log", params.toString());
		Job.animate();
		
		int count       = 0;
		int mxdp        = 20;
		int nkx2m       = 100*100;
		int nky2m       = nkx2m;
		double sfRP[][] = new double[(int)(Math.sqrt(nkx2m))+1][(int)(Math.sqrt(nky2m))+1];
		double sfIP[][] = new double[(int)(Math.sqrt(nkx2m))+1][(int)(Math.sqrt(nky2m))+1];
		double foo[]    = new double[2];
		
		while(count < mxdp){
			now = model.stepWT();
			Job.animate();
			if(now - then >= tau){
				then = now;
				for(int nx = 0; nx*nx <= nkx2m ; nx++){
					for(int ny = 0; ny*ny+nx*nx <= nky2m ; ny++){
						foo = model.structureFactor(nx, ny);
						sfRP[nx][ny] = foo[0];
						sfIP[nx][ny] = foo[1];
					}
					params.set("SF",Ef.format(100.*nx*nx/(double)(nkx2m)));
					Job.animate();
				}
				printSF(sfRP, sfIP, count, (int)(Math.sqrt(nkx2m))+1, (int)(Math.sqrt(nky2m))+1);
				count++;
			}
		}

		return;
	}
	
	private void printSF(double[][] sfRP, double[][] sfIP, int t, int xm, int ym){
		try{
			File file = new File(fout);
			PrintWriter pw = new PrintWriter(new FileWriter(file, true), true);
			if(t == 0)
				pw.println("t \t nx \t ny \t Re{S(k)} \t Im{S(k)}");
			for(int nx = 0; nx*nx <= xm ; nx++){
				for(int ny = 0; ny*ny <= ym ; ny++){
					if(sfRP[nx][ny] == 0)
						continue;
					pw.print(t +"\t");
					pw.print(nx+"\t");
					pw.print(ny+"\t");
					pw.println(sfRP[nx][ny]+"\t");
					pw.println(sfIP[nx][ny]);
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

package chris.MD.TwoD.apps;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.text.DecimalFormat;
import java.util.HashMap;
import java.util.Iterator;
import java.util.SortedSet;
import java.util.TreeSet;

import chris.MD.TwoD.LennardJones;
import chris.util.PrintUtil;
import scikit.jobs.Control;
import scikit.jobs.Job;
import scikit.jobs.Simulation;
import scikit.jobs.params.ChoiceValue;
import scikit.jobs.params.DirectoryValue;

public class LJradStrucfacApp extends Simulation{
	
	public LennardJones model;
	static DecimalFormat tf = new DecimalFormat("######.00");
	static DecimalFormat Ef = new DecimalFormat("0.###E0");
	private double now, then, tau, T;
	private String fout;
	
	
	public static void main(String[] args) {
		new Control(new LJradStrucfacApp(), "Lennard-Jones System");
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
		params.addm("d\u03C4",5);
		params.addm("T",1e-6);
		params.add("t");
		params.add("E");
		params.add("T_m");
		params.add("SF");

		
	}

	public void run() {
		
		double sf[][];
		int idx, countDG[];
		int countT                   = 0;
		int mxdp                     = 10;
		int nk2m                     = 1000*1000;
		SortedSet<Integer> ks2       = new TreeSet<Integer>();  
	    HashMap<Integer, Integer> hm = new HashMap<Integer, Integer>();
	    
		for(int jj = 0 ; (jj*jj) <= nk2m ; jj++){
			for(int kk = 0 ; (jj*jj + kk*kk) <= nk2m ; kk++){
				ks2.add(jj*jj+kk*kk);
			}
		}
		Iterator<Integer> it = ks2.iterator();
		while(it.hasNext()){
			hm.put(it.next(), countT++);
		}
		it     = null;
		countT = 0;
		
		sf      = new double[ks2.size()][mxdp];
		countDG = new int[ks2.size()];       

		
		model = new LennardJones(params);
		tau   = params.fget("d\u03C4");
		T     = params.fget("T");
		then  = 0;
		fout  = params.sget("Data Directory") + File.separator + params.sget("Data File") + ".txt";
		PrintUtil.printlnToFile(params.sget("Data Directory") + File.separator + "Params_" + params.sget("Data File") + ".log", params.toString());
		Job.animate();

		while(countT < mxdp){
			now = model.stepWT();
			Job.animate();
			if(now - then >= tau){
				then = now;
				for(int nx = 0; (nx*nx) <= nk2m ; nx++){
					for(int ny = 0; (idx = ny*ny+nx*nx) <= nk2m ; ny++){
						sf[hm.get(idx).intValue()][countT] = model.structureFactor(nx, ny);
						countDG[hm.get(idx).intValue()]++;
					}
					params.set("SF",Ef.format(100.*nx*nx/(double)(nk2m)));
					Job.animate();
				}
				for(int jj = 0 ; jj < countDG.length ; jj++){
					sf[jj][countT] /= countDG[jj];
					countDG[jj]     = 0;
				}
				countT++;
			}
		}

		printSF(sf, mxdp, ks2);
		return;
	}
	
	private void printSF(double[][] sf, int lngth, SortedSet<Integer> ss){
		
		Iterator<Integer> it;
		int count;

		try{
			File file = new File(fout);
			PrintWriter pw = new PrintWriter(new FileWriter(file, true), true);
			pw.println("t \t k^2 S(k)");
			for (int jj = 0 ; jj < lngth ; jj++){
				it = ss.iterator();
				count = 0;
				while(it.hasNext()){
					pw.print(jj+"\t");
					pw.print(it.next()+"\t");
					pw.println(sf[count++][jj]);
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

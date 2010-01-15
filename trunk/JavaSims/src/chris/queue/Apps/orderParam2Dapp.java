package chris.queue.Apps;

import java.io.File;
import java.text.DecimalFormat;

import scikit.dataset.Histogram;
import scikit.jobs.Control;
import scikit.jobs.Job;
import scikit.jobs.Simulation;
import scikit.jobs.params.DirectoryValue;
import scikit.jobs.params.DoubleValue;
import chris.queue.finiteRouter2D;
import chris.util.PrintUtil;

public class orderParam2Dapp extends Simulation{

	private DecimalFormat fmt = new DecimalFormat("0000");
	private int L, N, Mx, cycle, nOFt[];
	private String pth;
	private finiteRouter2D model;
	private Histogram h;

	public static void main(String[] args) {
		new Control(new orderParam2Dapp(), "2D Router Network");
	}
	
	public void load(Control c) {
		
		params.add("Data Directory",new DirectoryValue("/Users/cserino/Desktop/"));
		params.add("Data File", "default");
		params.add("M",100);
		params.add("L",16);
		params.add("l",5);
		params.add("\u03BB",new DoubleValue(0.05,0,1));
		params.add("seed",0);
		params.add("cycle");
		params.set("cycle",0);
		params.add("messages");
		params.set("messages",0);
		params.add("t");
		params.set("t",0);

		return;
	}

	public void run() {
		
		int tmx = (int)(1e5);
		nOFt    = new int[tmx+1];
		L       = params.iget("L");
		N       = L*L;
		Mx      = N*params.iget("M");
		pth     = params.sget("Data Directory");
		h       = new Histogram(1);

		for (cycle = 1 ; cycle < 1e5 ; cycle++){
			model = new finiteRouter2D(params); 
			params.set("cycle",cycle);
			Job.animate();
			while(model.step(true) < Mx && model.getT() < tmx){
				nOFt[model.getT()] += model.getNmsg();
			if(model.getT() % 1e3 == 0)
					Job.animate();
			}
			h.accum(model.getT());
			params.set("seed", params.iget("seed")+1);
			if (cycle % 100 == 0)
				save();
		}
	}
	
	public void save(){
		
		PrintUtil.printHistToFile(pth+File.separator+"hist"+fmt.format(cycle/100)+".txt", h);
		PrintUtil.printVectorToFile(pth+File.separator+"n"+fmt.format(cycle/100)+".txt", nOFt);
	}
	
	public void animate() {

		params.set("t", model.getT());
		params.set("messages",model.getNmsg());
		return;

	}

	public void clear() {
		
		return;
	}	
}

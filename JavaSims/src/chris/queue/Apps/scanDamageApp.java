package chris.queue.Apps;

import java.io.File;
import java.text.DecimalFormat;

import scikit.dataset.Histogram;
import scikit.jobs.Control;
import scikit.jobs.Job;
import scikit.jobs.Simulation;
import scikit.jobs.params.DirectoryValue;
import chris.queue.damagedRouter2D;
import chris.util.PrintUtil;

public class scanDamageApp extends Simulation{

	private String pth, fout;
	private damagedRouter2D model;
	private Histogram hN, ht;
	public static final DecimalFormat ifmt = new DecimalFormat("000");

	public static void main(String[] args) {
		new Control(new scanDamageApp(), "2D Router Network");
	}
	
	public void animate() {
		
		params.set("t", model.getT());
		params.set("messages",model.getNmsg());
		return;
	}

	public void clear() {
		
		return;
	}

	public void load(Control c) {
		
		params.add("Data Directory",new DirectoryValue("/Users/cserino/Desktop/"));
		params.add("Data File", "default");
		params.add("L");
		params.add("l");
		params.add("\u03BB");
		params.add("p");
		params.add("seed",0);
		params.add("messages");
		params.add("t");
		params.add("cycle");
		params.set("t",0);
		params.set("messages",0);
		params.set("L",32);


	}

	public void run() {
		pth   = params.sget("Data Directory");
		fout  = params.sget("Data File");
		
		double[] pv = new double[]{0., 0.01, 0.05, 0.1, 0.25}; 
		for(int ll = 5 ; ll < 25 ; ll += 5){
			params.set("l",ll);
			for(int aa = 0 ; aa < 20 ; aa++){
				params.set("\u03BB",1./(ll*(19.-aa)));
				for(int pp = 0 ; pp < 5 ; pp++){
					params.set("p",pv[pp]);
					model = new damagedRouter2D(params); 
					params.set("seed", params.iget("seed")+1);
					
					int count = 0;
					ht = new Histogram(1.);
					// equilibrate
					while(count < 1e4){
						model.step(1,false,ht);
						if((count++ % 500) == 0) 
							Job.animate();
					}		

					hN = new Histogram(1.);
					ht = new Histogram(1.);
					count = 0;
					// simulate model for data
					while(count < 1e6){
						model.step(1,false,ht);
						hN.accum(model.getNmsg());
						if((count++ % 500) == 0) 
							Job.animate();
					}
					printData(ll,aa,pp);
				}
			}
		}
	}

	private void printData(int x, int y, int z){
		String f1 = pth + File.separator + fout+"_nhist" +ifmt.format(x) +"_"+ifmt.format(y) +"_"+ifmt.format(z) +".txt";
		String f2 = pth + File.separator + fout+"_thist" +ifmt.format(x) +"_"+ifmt.format(y) +"_"+ifmt.format(z) +".txt";

		PrintUtil.printHistToFile(f1, ht);
		PrintUtil.printHistToFile(f2, hN);					
		return;
	}
}

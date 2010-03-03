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
	private static final int teql = 50000;
	private static final int tmax = 100000;
//	private static final int teql = 500;
//	private static final int tmax = 500;
	public static final DecimalFormat ifmt = new DecimalFormat("00");

	public static void main(String[] args) {
		new Control(new scanDamageApp(), "2D Router Network");
	}
	
	public void animate() {
		
		return;
	}

	public void clear() {
		
		return;
	}

	public void load(Control c) {
		
		params.add("Data Directory",new DirectoryValue("/Users/cserino/Desktop/"));
		params.add("Data File", "default");
		params.add("seed",0);
		params.add("L", 32);
		params.add("l", 5);
		params.add("\u03BB",0.19);
		params.add("p");
		params.add("messages");
		params.add("t");
		params.set("t",0);
		params.set("messages",0);


	}

	public void run() {

		double pp = 0;
		Histogram hN, ht;

		pth   = params.sget("Data Directory");
		fout  = params.sget("Data File");

		while(pp < 1){
			params.set("p",pp);
			params.set("seed", params.iget("seed")+1);
			model = new damagedRouter2D(params); 

			int count = 0;
			ht        = new Histogram(1.);
			// equilibrate
			while(count < teql){
				model.step(1,false,ht);
				if((count++ % 500) == 0){
					params.set("t", model.getT()-teql);
					params.set("messages",model.getNmsg());
					Job.animate();
				}
			}

			count = 0;
			int[] ts  = new int[tmax];
			hN = new Histogram(1.);
			ht = new Histogram(1.);
			// simulate model for data
			while(count < tmax){
				model.step(1,false,ht);
				hN.accum(model.getNmsg());
				ts[count] = model.getNmsg();
				if((count++ % 500) == 0){ 
					params.set("t", model.getT());
					params.set("messages",model.getNmsg());
					Job.animate();
				}
			}
			printData(pp, hN, ht, ts);

			if(pp < 0.1){
				pp = Math.round(100*(pp+0.01))/100.; 
			}
			else if(pp < 0.2){
				pp = Math.round(100*(pp+0.1))/100.; 
			}
			else{
				pp = Math.round(100*(pp+0.2))/100.; 
			}
		}
	}


	private void printData(double pp, Histogram h1, Histogram h2, int[] v){
		String f1 = pth + File.separator + fout+"_nhist_" +ifmt.format(100*pp)+".txt";
		String f2 = pth + File.separator + fout+"_thist_" +ifmt.format(100*pp)+".txt";
		String f3 = pth + File.separator + fout+"_tsers_" +ifmt.format(100*pp)+".txt";

		PrintUtil.printHistToFile(f1, h1);
		PrintUtil.printHistToFile(f2, h2);	
		PrintUtil.printVectorToFile(f3,v);
		return;
	}
}

package chris.queue.Apps;

import java.io.File;
import java.text.DecimalFormat;

import scikit.dataset.Histogram;
import scikit.jobs.Control;
import scikit.jobs.Job;
import scikit.jobs.Simulation;
import scikit.jobs.params.DirectoryValue;
import scikit.jobs.params.DoubleValue;
import chris.queue.damagedRouter2D;
import chris.util.MathUtil;
import chris.util.PrintUtil;

public class damageScaling2Dapp extends Simulation{
	
	public static final DecimalFormat ifmt = new DecimalFormat("00");
	private static final int teql = 100000;
	private static final int tmax = 500000;
	private String pth, fout;
	private damagedRouter2D model;

	public static void main(String[] args) {
		new Control(new damageScaling2Dapp(), "2D Router Network");
	}
	
	public void load(Control c) {
		
		params.add("Data Directory",new DirectoryValue("/Users/cserino/Desktop/"));
		params.add("Data File", "default");
		params.add("seed",0);
		params.add("L",32);
		params.add("l",5);
		params.add("p",new DoubleValue(0.,0,1));
		params.add("\u03BB");
		params.add("messages");
		params.add("t");
		params.add("cycle");
	}

	public void run() {
		
		pth   = params.sget("Data Directory");
		fout  = pth + File.separator + params.sget("Data File");
		
		double[] lv = MathUtil.logspace(Math.log10(0.19), Math.log10(0.199),25);
		PrintUtil.printVectorToFile(fout+"_jj2lambda.txt", lv);
		PrintUtil.printlnToFile(pth + File.separator + "Params_" + params.sget("Data File") + ".log", params.toString());

		for(int jj = 2 ; jj < lv.length; jj++){
			params.set("cycle",jj+1);
			params.set("seed", params.iget("seed") + 1);
			params.set("\u03BB",lv[jj]);
			model = new damagedRouter2D(params); 

			Histogram hnmsg = new Histogram(1.);
			Histogram htlat = new Histogram(1.);
			int count = 0;
			while(count < teql){
				model.step(1,false,htlat);
				if(count++ % 1000 == 0){
					params.set("t", model.getT()-1-teql);
					params.set("messages",model.getNmsg());
					Job.animate();
				}
			}

			count = 0;
			hnmsg = new Histogram(1.);
			htlat = new Histogram(1.);
			while(count < tmax){
				model.step(1,false,htlat);
				hnmsg.accum(model.getNmsg());
				if(count++ % 1000 == 0){
					params.set("t", model.getT()-1);
					params.set("messages",model.getNmsg());
					Job.animate();
				}
			}
			printData(jj, hnmsg, htlat);
			model = null;
		}
		
		params.set("cycle","done");
		Job.animate();
		return;
	}
	
	private void printData(int index, Histogram h1, Histogram h2){
		String f1 = fout+"_nhist" +ifmt.format(index)+".txt";
		String f2 = fout+"_thist" +ifmt.format(index)+".txt";

		PrintUtil.printHistToFile(f1, h1);
		PrintUtil.printHistToFile(f2, h2);	
		return;
	}
	
	
	public void animate() {}

	public void clear() {}
}

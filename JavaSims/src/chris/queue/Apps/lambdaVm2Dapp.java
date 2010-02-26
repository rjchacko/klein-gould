package chris.queue.Apps;

import java.io.File;

import scikit.dataset.Histogram;
import scikit.jobs.Control;
import scikit.jobs.Job;
import scikit.jobs.Simulation;
import scikit.jobs.params.DirectoryValue;
import chris.queue.finiteRouter2D;
import chris.util.PrintUtil;

public class lambdaVm2Dapp extends Simulation{

	private int L, N, M, Mx, tw;
	private double lambda;
	private String pth;
	private finiteRouter2D model;

	public static void main(String[] args) {
		new Control(new lambdaVm2Dapp(), "2D Router Network");
	}
	
	public void load(Control c) {
		
		params.add("Data Directory",new DirectoryValue("/Users/cserino/Desktop/"));
		params.add("Data File", "default");
		params.add("seed",0);
		params.add("L",32);
		params.add("l",5);
		params.add("M",10);
		params.add("t_wait",(int)(1e6));
		params.add("messages");
		params.set("messages",0);
		params.add("t");
		params.set("t",0);
		params.add("\u03BB");
	}

	public void run() {
	
		int count;
		Histogram foo = new Histogram(1000.);

		
		L   = params.iget("L");
		N   = L*L;
		pth = params.sget("Data Directory") + File.separator + params.sget("Data File") + ".txt";
		tw  = params.iget("t_wait");
		M   = params.iget("M");
		PrintUtil.printlnToFile(pth,"M","lambda");

		while(true){ // loop on m
		
			params.set("M",M);
			Mx = M*N;	
			//lambda = 0.1999;
			lambda = 0.12;
			
			while(true){ // loop on lambda
				params.set("\u03BB",lambda);
				model = null; // avoid heap space issues
				model = new finiteRouter2D(params); 
				params.set("seed",params.iget("seed")+1);
				Job.animate();
				count = 0;
				while(model.step(1,false,foo) < Mx && count++ < tw)
					maybeAnimate(count);

				params.set("t","Done");
				if(count >= tw)
					break;
				
				if(lambda == 0.19999){
					lambda = 0.1999;
				}
				else if(lambda == 0.1999){
					lambda = 0.199;
				}
				else{
					lambda -= 0.001;
					lambda = Math.round(1000*lambda)/1000.; // avoid annoying 0.500000000001
				}
			}
			// found lambda s.t. system did not nucleate in a given time
			PrintUtil.printlnToFile(pth,M,lambda);
			int p = (int)(Math.floor(Math.log10(M)));
			M     += (int)(Math.pow(10, p));
		}
	}
	
	private void maybeAnimate(int now){
		
		if(now%1000 == 0)
			Job.animate();
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

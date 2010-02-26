package chris.queue.Apps;

import java.io.File;

import scikit.dataset.Histogram;
import scikit.jobs.Control;
import scikit.jobs.Job;
import scikit.jobs.Simulation;
import scikit.jobs.params.DirectoryValue;
import chris.queue.finiteRouter2D;
import chris.util.PrintUtil;

public class susceptibility2Dapp extends Simulation{

	private int L, N, Mx;
	private String pth;
	private finiteRouter2D model;
	private Histogram hN, ht;


	public static void main(String[] args) {
		new Control(new susceptibility2Dapp(), "2D Router Network");
	}
	
	public void load(Control c) {
		
		params.add("Data Directory",new DirectoryValue("/Users/cserino/Desktop/"));
		params.add("Data File", "default");
		params.add("M");
		params.add("L");
		params.add("l");
		params.add("\u03BB");
		params.add("seed",0);
		params.add("messages");
		params.add("t");
		params.add("cycle");
		params.set("t",0);
		params.set("messages",0);
		params.set("M",100);
		params.set("L",32);
		params.set("l",5);

	}

	public void run() {
		
		L     = params.iget("L");
		N     = L*L;
		Mx    = N*params.iget("M");
		pth   = params.sget("Data Directory");
		double[] lambdav = new double[]{0.000100000 ,0.000106886 ,0.000114247 ,0.000122115 ,0.000130524 ,0.000139512 ,0.000149120 ,0.000159389 ,0.000170365 ,0.000182097 ,0.000194637 ,0.000208040 ,0.000222367 ,0.000237680 ,0.000254047 ,0.000271542 ,0.000290241 ,0.000310229 ,0.000331592 ,0.000354427 ,0.000378834 ,0.000404922 ,0.000432807 ,0.000462611 ,0.000494469 ,0.000528520 ,0.000564916 ,0.000603818 ,0.000645400 ,0.000689844 ,0.000737350 ,0.000788127 ,0.000842400 ,0.000900411 ,0.000962417 ,0.001028693 ,0.001099533 ,0.001175251 ,0.001256184 ,0.001342690 ,0.001435153 ,0.001533983 ,0.001639619 ,0.001752530 ,0.001873216 ,0.002002214 ,0.002140094 ,0.002287469 ,0.002444994 ,0.002613366 ,0.002793333 ,0.002985693 ,0.003191299 ,0.003411065 ,0.003645965 ,0.003897041 ,0.004165406 ,0.004452253 ,0.004758853 ,0.005086567 ,0.005436848 ,0.005811251 ,0.006211437 ,0.006639182 ,0.007096382 ,0.007585068 ,0.008107406 ,0.008665714 ,0.009262470 ,0.009900321 ,0.010582097 ,0.011310822 ,0.012089731 ,0.012922278 ,0.013812158 ,0.014763318 ,0.015779980 ,0.016866652 ,0.018028158 ,0.019269649 ,0.020596634 ,0.022015001 ,0.023531042 ,0.025151484 ,0.026883516 ,0.028734822 ,0.030713617 ,0.032828680 ,0.035089395 ,0.037505791 ,0.040088590 ,0.042849251 ,0.045800022 ,0.048953995 ,0.052325164 ,0.055928484 ,0.059779944 ,0.063896630 ,0.068296808 ,0.073000000};
		boolean nuc = false;
		PrintUtil.printlnToFile(pth+File.separator+"Params_"+params.sget("Data File")+".txt", params.toString());
		PrintUtil.printVectorToFile(pth+File.separator+"jj2lambda.txt", lambdav);
		for(int jj = 0 ; jj < lambdav.length ; jj++){
			params.set("cycle",jj);
			params.set("\u03BB",lambdav[jj]);
			model = new finiteRouter2D(params); 
			params.set("seed", params.iget("seed")+1);
			int count = 0;
			// equilibrate model
			hN = new Histogram(1.);
			ht = new Histogram(1.);
			while(model.step(1,false,ht) < Mx && count < 1e4){
				hN.accum(model.getNmsg());
				if((count++ % 500) == 0) 
					Job.animate();
			}		
			count = 0;
			hN = new Histogram(1.);
			ht = new Histogram(1.);
			// simulate model for data
			while(model.step(1,false,ht) < Mx && count < 1e6){
				hN.accum(model.getNmsg());
				if((count++ % 500) == 0) 
					Job.animate();
			}
			nuc = (model.getNmsg() >= Mx); 
			// write data here
			if(nuc){
				PrintUtil.printHistToFile(pth+File.separator+"tfHist_"+params.sget("Data File")+jj+"XXX.txt", ht);
				PrintUtil.printHistToFile(pth+File.separator+"nHist_"+params.sget("Data File")+jj+"XXX.txt", hN);
			}
			else{
				PrintUtil.printHistToFile(pth+File.separator+"tfHist_"+params.sget("Data File")+jj+".txt", ht);
				PrintUtil.printHistToFile(pth+File.separator+"nHist_"+params.sget("Data File")+jj+".txt", hN);				
			}
			nuc   = false;
			model = null;
		}


		params.set("t","Done");
		Job.signalStop();
		Job.animate();
		return;
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

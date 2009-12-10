package chris.queue.Apps;

import scikit.jobs.Control;
import scikit.jobs.Job;
import scikit.jobs.Simulation;
import scikit.jobs.params.DirectoryValue;
import scikit.jobs.params.DoubleValue;
import chris.queue.router1D;

public class routerData1Dapp extends Simulation{

	int N, tmax;
	router1D model;

	public static void main(String[] args) {
		new Control(new routerData1Dapp(), "1D Router Network");
	}
	
	public void load(Control c) {
		
		params.add("Data Directory",new DirectoryValue("/Users/cserino/Desktop/"));
		params.add("Data File", "default");
		params.add("N",50);
		params.add("l",10);
		params.add("\u03BB",new DoubleValue(0.05,0,1));
		params.add("seed",0);
		params.add("messages");
		params.set("messages",0);
		params.add("t_max",(int)(1e8));
		params.add("t");
		params.set("t",0);
	}

	public void run() {
		
		int count = 0;
		
		N     = params.iget("N");
		model = new router1D(params); 
		tmax  = params.iget("t_max");
			
		while(count++ < tmax){
			model.step();
			if(count % 10000 == 0)
				Job.animate();
		}
		
		Job.signalStop();
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

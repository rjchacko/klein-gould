package chris.queue.Apps;

import scikit.jobs.Control;
import scikit.jobs.Job;
import scikit.jobs.Simulation;
import scikit.jobs.params.DirectoryValue;
import scikit.jobs.params.DoubleValue;
import chris.queue.finiteRouter2D;

public class finiteRouterData2Dapp extends Simulation{

	private int L, N, Mx;
	private finiteRouter2D model;

	public static void main(String[] args) {
		new Control(new finiteRouterData2Dapp(), "2D Router Network");
	}
	
	public void load(Control c) {
		
		params.add("Data Directory",new DirectoryValue("/Users/cserino/Desktop/"));
		params.add("Data File", "default");
		params.add("M",10);
		params.add("L",32);
		params.add("l",5);
		params.add("\u03BB",new DoubleValue(0.05,0,1));
		params.add("seed",0);
		params.add("cycle");
		params.set("cycle",0);
		params.add("messages");
		params.set("messages",0);
		params.add("t");
		params.set("t",0);

	}

	public void run() {
		
		L     = params.iget("L");
		N     = L*L;
		Mx    = N*params.iget("M");
		model = new finiteRouter2D(params); 
		
		
		for(int jj = 0 ; jj < 100 ; jj++){
			model.step(false);
		}
		
		while(model.step(true) < Mx)
			if(model.getT() % 1000 == 0) 
				Job.animate();
		model.writeData(model.getT());
		model.writePPdata(model.getT(), params.iget("cycle"));
		
		params.set("t","Done");
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

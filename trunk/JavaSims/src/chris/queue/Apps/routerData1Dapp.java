package chris.queue.Apps;

import scikit.graphics.dim2.Grid;
import scikit.jobs.Control;
import scikit.jobs.Job;
import scikit.jobs.Simulation;
import scikit.jobs.params.DirectoryValue;
import chris.queue.router1D;

public class routerData1Dapp extends Simulation{

	Grid grid = new Grid("Buffers");
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
		params.add("\u03BB",0.05);
		params.add("seed",0);
		params.add("messages");
		params.set("messages",0);
		params.add("t_max",1e7);
		params.add("t");
		params.set("t",0);

		c.frame(grid);
	}

	public void run() {
		
		N     = params.iget("N");
		model = new router1D(params); 
		int count = 0;
		
		while(count++ < tmax){
			model.step();
			if(count % 5000 == 0)
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
		
		grid.clear();
		return;
	}	
}

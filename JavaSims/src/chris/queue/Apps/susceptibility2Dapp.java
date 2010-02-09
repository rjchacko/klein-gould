package chris.queue.Apps;

import java.io.File;

import scikit.dataset.Histogram;
import scikit.jobs.Control;
import scikit.jobs.Job;
import scikit.jobs.Simulation;
import scikit.jobs.params.DirectoryValue;
import scikit.jobs.params.DoubleValue;
import chris.queue.finiteRouter2D;
import chris.util.PrintUtil;

public class susceptibility2Dapp extends Simulation{

	private int L, N, Mx;
	private String pth;
	private finiteRouter2D model;


	public static void main(String[] args) {
		new Control(new susceptibility2Dapp(), "2D Router Network");
	}
	
	public void load(Control c) {
		
		params.add("Data Directory",new DirectoryValue("/Users/cserino/Desktop/"));
		params.add("Data File", "default");
		params.add("M",10);
		params.add("L",32);
		params.add("l",5);
		params.add("\u03BB",new DoubleValue(0.1,0,1));
		params.add("seed",0);
		params.add("messages");
		params.set("messages",0);
		params.add("t");
		params.set("t",0);
		params.add("cycle");

	}

	public void run() {
		
		L     = params.iget("L");
		N     = L*L;
		Mx    = N*params.iget("M");
		pth   = params.sget("Data Directory");
		
		double lambda = params.fget("\u03BB");
		while(lambda < 0.145){
			Histogram h = new Histogram(1);
			for(int jj = 0 ; jj < 100 ; jj++){
				params.set("cycle",jj);
				model = new finiteRouter2D(params); 
				params.set("seed", params.iget("seed")+1);
				int count = 0;
				while(model.step(false) < Mx && count < 2e3)
					if((count++ % 500) == 0) 
						Job.animate();
				h.accum(model.getT());
				model = null;
			}
			lambda += 0.005;
			params.set("\u03BB",lambda);
			PrintUtil.printHistToFile(pth+File.separator+"tfHist"+finiteRouter2D.cmt.format(lambda*1000)+".txt", h);
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

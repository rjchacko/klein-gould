package chris.RandomWalker;

import java.io.File;

import scikit.jobs.Control;
import scikit.jobs.Job;
import scikit.jobs.Simulation;
import scikit.jobs.params.DirectoryValue;
import scikit.jobs.params.DoubleValue;
import scikit.jobs.params.StringValue;

public class PearsonApp extends Simulation{

	private int d, nw, ns;
	private double lambda;
	private ShrinkingPearson model;
	
	public static void main(String[] args) {
		new Control(new PearsonApp(), "Parameters");
	}
	
	public void load(Control c) {
	
		params.add("Data Directory",new DirectoryValue("/Users/cserino/Desktop"));
		params.add("File Name", new StringValue("default"));
		params.add("Random Seed",(int) 0);
		params.add("Dimension",(int) 2);
		params.add("Resolution", (double) 0.01);
		params.add("Walkers", (int) 1e8);
		params.add("Steps", (int) 10);
		params.add("lambda", new DoubleValue(0.5, 0.5, 1.));
		params.add("Status");
	}

	public void run() {
		
		nw     = params.iget("Walkers");
		ns     = params.iget("Steps");
		d      = params.iget("Dimension");
		lambda = params.fget("lambda");
		model  = new ShrinkingPearson(params);
		
		for (int jj = 0 ; jj < nw ; jj++){
			model.nextWalk(lambda, ns, d);
			if(jj % 500 == 0){
				params.set("Status",jj);
				Job.animate();
			}
		}
		
		model.printPDF(params.sget("Data Directory")+File.separator+params.sget("File Name")+".txt", ns, nw, lambda);
		
		params.set("Status","Done");
		Job.animate();
	}
	
	public void animate() {
		
		return;
	}

	public void clear() {
	
		return;
	}

	
}

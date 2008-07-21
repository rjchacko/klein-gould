package chris.misc;

import chris.util.PrintUtil;
import scikit.jobs.Control;
import scikit.jobs.Job;
import scikit.jobs.Simulation;

public class WriteCourse extends Simulation{

	
	
	public static void main(String[] args) {
		new Control(new WriteCourse(), "Golf Score");
	}
	
	
	public void animate() {
		
	}

	public void clear() {
		params.set("Distance", (int) 0);
		params.set("Par", (int) 0);
		params.set("HCP", (int) 0);
	}

	public void load(Control c) {
		
		params.add("Hole");
		params.set("Hole","Press 'Step' to Begin");
		params.addm("Distance", (int) 0);
		params.addm("Par", (int) 0);
		params.addm("HCP", (int) 0);
	}

	public void run() {

		int dist, par, hcp;
		int hole = 1;
		
		params.set("Hole",hole++);

		Job.animate();

		
		dist = params.iget("Distance");
		par  = params.iget("Par");
		hcp  = params.iget("HCP");
		
		PrintUtil.printlnToFile("/home/cserino/Desktop/test2.txt","\\addscore{"+(hole-1)+"}{"+dist+"}{"+par+"}{"+hcp+"}");	

		clear();
		params.set("Hole",hole++);	
		Job.animate();

		while(hole <= 19){
			
			dist = params.iget("Distance");
			par  = params.iget("Par");
			hcp  = params.iget("HCP");

			PrintUtil.printlnToFile("/home/cserino/Desktop/test2.txt","\\addscore{"+(hole-1)+"}{"+dist+"}{"+par+"}{"+hcp+"}");	
			
			clear();
			params.set("Hole",hole++);	
			if(hole > 19) params.set("Hole","Done");	
			Job.animate();
			
		}
		
		return;
		
	}

}

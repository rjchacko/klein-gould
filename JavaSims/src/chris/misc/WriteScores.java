package chris.misc;

import chris.util.PrintUtil;
import scikit.jobs.Control;
import scikit.jobs.Job;
import scikit.jobs.Simulation;

public class WriteScores extends Simulation{

	
	
	public static void main(String[] args) {
		new Control(new WriteScores(), "Golf Score");
	}
	
	
	public void animate() {
		
	}

	public void clear() {
		
		params.set("Score", (int) 0);
		params.set("+/-", (int) 0);
	}

	public void load(Control c) {
		
		params.add("Hole");
		params.set("Hole","Press 'Step' to Begin");
		params.addm("Score", (int) 0);
		params.addm("+/-", (int) 0);
	}

	public void run() {

		int score, pm;
		int hole = 1;
		

		while(hole < 19){
			
			params.set("Hole",hole++);
			Job.animate();
			
			score = params.iget("Score");
			pm    = params.iget("+/-");
			
			PrintUtil.printlnToFile("/Users/cserino/Desktop/test2.txt","\\addscore{"+(hole-1)+"}{"+score+"}{"+pm+"}");	

			clear();
			params.set("Hole",hole);	
			if(hole > 18){
				params.set("Hole","Done");	
				Job.animate();
				return;
			}
			
		}
		
		return;
		
	}

}

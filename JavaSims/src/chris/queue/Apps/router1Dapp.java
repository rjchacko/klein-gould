package chris.queue.Apps;

import scikit.graphics.dim2.Grid;
import scikit.jobs.Control;
import scikit.jobs.Job;
import scikit.jobs.Simulation;
import chris.queue.router1D;

public class router1Dapp extends Simulation{

	Grid grid = new Grid("Buffers");
	int N;
	router1D model;
//	ColorGradient cg;

	public static void main(String[] args) {
		new Control(new router1Dapp(), "1D Router Network");
	}
	
	public void load(Control c) {
		
		params.add("N",20);
		params.add("l",5);
		params.add("\u03BB",0.01);
		params.add("seed",0);
		params.add("messages");
		params.set("messages",0);
		params.add("t");
		params.set("t",0);

		c.frame(grid);
	}

	public void run() {
		
		N     = params.iget("N");
		model = new router1D(params); 
		
		// H E R E !!!
	//	setupcolors();
//		cg    = new ColorGradient();
//		int v = new int[]
//		grid.setColors(cg.getColor(v,0,5*lambda*l)
				
		while(true){
			model.step();
			Job.animate();
		}
		
	}

	public void animate() {
		
//		ColorPalette palette = new ColorPalette();
//		palette.setColor(0,Color.WHITE);
//		palette.setColor(1,Color.BLACK);
//		palette.setColor(2,Color.RED);
//		palette.setColor(3,Color.BLUE);
//		grid.setColors(palette);
		
		params.set("t", model.getT());
		params.set("messages",model.getNmsg());
		grid.registerData(N,1,model.getDensity());
		return;
	}

	public void clear() {
		
		grid.clear();
		return;
	}	
}

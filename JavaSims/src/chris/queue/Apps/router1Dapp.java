package chris.queue.Apps;

import scikit.graphics.ColorGradient;
import scikit.graphics.ColorPalette;
import scikit.graphics.dim2.Grid;
import scikit.jobs.Control;
import scikit.jobs.Job;
import scikit.jobs.Simulation;
import scikit.jobs.params.DirectoryValue;
import scikit.jobs.params.DoubleValue;
import chris.queue.router1D;

public class router1Dapp extends Simulation{

	Grid grid = new Grid("Buffers");
	int N;
	router1D model;
	private static final boolean draw = true;

	public static void main(String[] args) {
		new Control(new router1Dapp(), "1D Router Network");
	}
	
	public void load(Control c) {
		
		params.add("Data Directory",new DirectoryValue("/Users/cserino/Desktop/"));
		params.add("Data File", "default");
		params.add("N",20);
		params.add("l",5);
		params.add("\u03BB",new DoubleValue(0.05,0,1));
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

		setupcolors(params.fget("\u03BB"),params.iget("l"),N);
				
		while(true){
			model.step(true);
			Job.animate();
		}
		
	}
	
	private void setupcolors(double lambda, int l, int L){

		int nmbar = (int)((lambda*lambda*(l-1)/(1-l*lambda) + l*lambda));
		if(nmbar == 0){
			nmbar = 5;
		}
		else{
			nmbar = 2*(nmbar+1)*N;
		}
		setupcolors(nmbar);
		return;
	}

	private void setupcolors(int mx){
		
		ColorGradient cg = new ColorGradient();
		ColorPalette cp  = new ColorPalette();

		for (int jj = 0 ; jj <= mx ; jj++){
			cp.setColor(mx-jj,cg.getColor(jj, 0, mx));
		}
		grid.setColors(cp);
		return;
	}
	
	public void animate() {
		
		params.set("t", model.getT());
		params.set("messages",model.getNmsg());
		if(draw) 
			grid.registerData(N,1,model.getDensity());
		return;
	}

	public void clear() {
		
		grid.clear();
		return;
	}	
}

package rachele.ising.dim2MC.apps;

import static scikit.util.Utilities.format;

import rachele.ising.dim2MC.IsingLRCircle;
import rachele.util.FourierTransformer;
import scikit.graphics.dim2.Grid;
import scikit.jobs.Control;
import scikit.jobs.Job;
import scikit.jobs.Simulation;
import scikit.jobs.params.ChoiceValue;

public class MCIsingCircleApp extends Simulation{
	
	Grid grid = new Grid("Long Range Ising Model");
	Grid sfGrid = new Grid("SF");
	IsingLRCircle sim;
	FourierTransformer fft;
	boolean measure=false;
	double [] ft;

	
	public static void main(String[] args) {
		new Control(new MCIsingCircleApp(), "Monte Carlo");
	}

	public void load(Control c) {
		c.frameTogether("data",grid, sfGrid);
		params.addm("Dynamics", new ChoiceValue("Ising Glauber","Kawasaki Glauber", "Kawasaki Metropolis",  "Ising Metropolis"));
		params.add("Random seed", 0);
		params.add("L", 1<<5);
		params.add("R", 10);//1<<6);
		params.add("Initial magnetization", 0.0);
		params.addm("T", 0.04);
		params.addm("J", -1.0);
		params.addm("h", 0.0);
		params.addm("dt", 1.0);//1/(double)(1<<4));
		params.add("time");
		params.add("magnetization");
	}
	

	public void animate() {
		grid.setScale(-1.0, 1.0);
		grid.registerData(sim.L, sim.L, sim.spins.spin);
		sfGrid.registerData(sim.L, sim.L, ft);
		params.set("time", format(sim.time()));
		sim.setParameters(params);
			
	}


	public void clear() {
	}


	public void run() {
		sim = new IsingLRCircle(params);
		sim.randomizeField(params.fget("Initial magnetization"));	
		fft = new FourierTransformer(sim.L);
		ft = new double [sim.L*sim.L];
		while(true){
			sim.step();
			ft = fft.calculate2DSF(sim.getField(), true, true);
			Job.animate();
		}
		
	}

}

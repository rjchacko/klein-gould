package rachele.ising.dim2MC.apps;

import static scikit.util.Utilities.format;

import java.awt.Color;

import rachele.ising.dim2MC.IsingLRCircle;
import rachele.util.FourierTransformer;
import scikit.dataset.Accumulator;
import scikit.graphics.dim2.Grid;
import scikit.graphics.dim2.Plot;
import scikit.jobs.Control;
import scikit.jobs.Job;
import scikit.jobs.Simulation;
import scikit.jobs.params.ChoiceValue;

public class MCIsingCircleApp extends Simulation{
	
	Grid grid = new Grid("Long Range Ising Model");
	IsingLRCircle sim;
	FourierTransformer fft;
	Plot sf = new Plot("SF");
	Plot sf2 = new Plot("SF");
	boolean measure=false;
	double [] ft;
	Accumulator sfh = new Accumulator();
	Accumulator sfv = new Accumulator();
	
	
	public static void main(String[] args) {
		new Control(new MCIsingCircleApp(), "Monte Carlo");
	}

	public void load(Control c) {
		c.frameTogether("data",grid, sf, sf2);
		params.addm("Dynamics", new ChoiceValue("Ising Glauber","Kawasaki Glauber", "Kawasaki Metropolis",  "Ising Metropolis"));
		params.add("Random seed", 0);
		params.add("L", 1<<7);
		params.add("R", 32);//1<<6);
		params.add("Initial magnetization", 0.0);
		params.addm("T", 0.04);
		params.addm("J", -1.0);
		params.addm("h", 0.0);
		params.addm("dt", 1.0);//1/(double)(1<<4));
		params.addm("k int", 3);
		params.add("time");
		params.add("magnetization");
		flags.add("Clear");
	}
	

	public void animate() {
		grid.setScale(-1.0, 1.0);
		grid.registerData(sim.L, sim.L, sim.spins.spin);
		sf.setAutoScale(true);
		sf2.setAutoScale(true);
		sf.registerPoints("sf h", sfh, Color.BLUE);
		sf2.registerPoints("sf v", sfv, Color.RED);
		params.set("time", format(sim.time()));
		sim.setParameters(params);
		if (flags.contains("Clear")){
			sfh.clear();
			sfv.clear();
			flags.clear();
		}
			
	}


	public void clear() {
		sfh.clear();
		sfv.clear();
	}


	public void run() {
		sim = new IsingLRCircle(params);
		sim.randomizeField(params.fget("Initial magnetization"));	
		fft = new FourierTransformer(sim.L);
		ft = new double [sim.L*sim.L];
		int kInt = params.iget("k int");
		System.out.println("k = " + 2*Math.PI*kInt*sim.R/sim.L + " spin o " + sim.spins.spin[0]);
		while(true){
			sim.step();
//			ft = fft.calculate2DSF(sim.getField(), false, false);
//			sfh.accum(sim.time(), ft[kInt]);
//			sfh.accum(sim.time(), ft[(-kInt+sim.L)%sim.L]);
//			sfv.accum(sim.time(), ft[kInt*sim.L]);
//			sfv.accum(sim.time(), ft[((-kInt+sim.L)%sim.L)*sim.L]);
			Job.animate();
		}
		
	}

}

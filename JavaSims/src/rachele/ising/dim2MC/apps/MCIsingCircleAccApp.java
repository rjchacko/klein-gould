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

public class MCIsingCircleAccApp extends Simulation{
	
	Grid grid = new Grid("Long Range Ising Model");
	Grid sfGrid = new Grid("SF");
	Plot sfs = new Plot("Sfs");
	IsingLRCircle sim;
	FourierTransformer fft;
	boolean measure=false;
	double [] ft;
	Accumulator sfh = new Accumulator();
	Accumulator sfv = new Accumulator();
	Accumulator sf0 = new Accumulator();
	
	
	public static void main(String[] args) {
		new Control(new MCIsingCircleAccApp(), "Monte Carlo");
	}

	public void load(Control c) {
		c.frameTogether("data",grid, sfGrid, sfs);
		params.addm("Dynamics", new ChoiceValue("Ising Glauber","Kawasaki Glauber", "Kawasaki Metropolis",  "Ising Metropolis"));
		params.add("Random seed", 0);
		params.add("L", 1<<5);
		params.add("R", 10);//1<<6);
		params.add("Initial magnetization", 0.0);
		params.addm("T", 0.04);
		params.addm("J", -1.0);
		params.addm("h", 0.0);
		params.addm("dt", 1.0);//1/(double)(1<<4));
		params.addm("kx int", 3);
		params.addm("ky int", 3);
		params.add("time");
		params.add("magnetization");
		flags.add("Clear");
	}

	public void animate() {
		grid.setScale(-1.0, 1.0);
		grid.registerData(sim.L, sim.L, sim.spins.spin);
		sfGrid.registerData(sim.L, sim.L, ft);
		params.set("time", format(sim.time()));
		sim.setParameters(params);
		sfs.registerPoints("Sf0", sf0, Color.BLACK);
		sfs.registerPoints("SfH", sfh, Color.BLUE);
		sfs.registerPoints("SfV", sfv, Color.RED);
		if (flags.contains("Clear")){
			sfh.clear();
			sfv.clear();
			sf0.clear();
			flags.clear();
		}
			
	}

	public void clear() {
		sfh.clear();
		sfv.clear();
		sf0.clear();
	}

	public void run() {
		sim = new IsingLRCircle(params);
		sim.randomizeField(params.fget("Initial magnetization"));	
		fft = new FourierTransformer(sim.L);
		ft = new double [sim.L*sim.L];
		int kxInt = params.iget("kx int");
		int kyInt = params.iget("ky int");
		System.out.println("kx = " + 2*Math.PI*kxInt*sim.R/sim.L);
		System.out.println("ky = " + 2*Math.PI*kyInt*sim.R/sim.L);
		while(true){
			sim.step();
			ft = fft.calculate2DSF(sim.getField(), false, false);
			sf0.accum(sim.time(), ft[0]);
			sfh.accum(sim.time(), ft[kxInt]);
			sfh.accum(sim.time(), ft[(sim.L-kxInt)%sim.L]);
			sfv.accum(sim.time(), ft[kyInt*sim.L]);
			sfv.accum(sim.time(), ft[((sim.L-kxInt)%sim.L)*sim.L]);
			Job.animate();
		}
	}

}
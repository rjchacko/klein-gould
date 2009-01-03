package kip.ising.dim2.apps;

import kip.ising.NewmanZiff;
import kip.util.Random;
import scikit.graphics.dim2.Grid;
import scikit.jobs.Control;
import scikit.jobs.Job;
import scikit.jobs.Simulation;
import scikit.util.Utilities;


public class Ising2DApp extends Simulation {
    Grid grid = new Grid("Ising spins");
    Grid perc = new Grid("Perc");
	Ising sim;
	double dt = 0.5;
	double[] clusters;
	
	public static void main(String[] args) {
		new Control(new Ising2DApp(), "Ising Model");
	}

	public void load(Control c) {
		c.frame(grid, perc);
		
		params.add("Seed", 0);
		params.add("L", 256);
		params.add("Ratio", 1.0);
		params.add("T", 0.0);
		params.addm("dt", 1.0);
		params.add("time");
		params.add("homology");
	}
	
	public void animate() {
		params.set("time", Utilities.format(sim.time));
		
		sim.T = params.fget("T");
		dt = params.fget("dt");
		grid.registerData(sim.L1, sim.L2, sim.spin);
		
		NewmanZiff nz = sim.findClusters(1);
		for (int i = 0; i < clusters.length; i++)
			clusters[i] = nz.clusterSize(i);
		perc.registerData(sim.L1, sim.L2, clusters);
		params.set("homology", (nz.wrap_horizontal ? "horiz ":"") + (nz.wrap_vertical ? "vert ":"") + (nz.wrap_diagonal ? "diag ":""));
	}
	
	public void clear() {
		grid.clear();
	}
	
	public void run() {
		int seed = params.iget("Seed");
		int L1 = params.iget("L");
		int L2 = (int) (L1 * params.fget("Ratio"));
		sim = new Ising(seed, L1, L2, params.fget("T"));	
		clusters = new double[L1*L2];
		
        while (true) {
        	sim.step(dt);            
        	Job.animate();
        }
	}
	

	class Ising {
		public int spin[];
		public int L1, L2;
		public int N;
		public double T;
		public Random random = new Random();
		public double time;
		
		public Ising(int seed, int _L1, int _L2, double _T) {
			random.setSeed(seed);
			L1 = _L1;
			L2 = _L2;
			T = _T;
			N = L1*L2;
			time = 0;
			
			spin = new int[N];
			randomize();
		}
		
		public void randomize() {
			for (int i = 0; i < N; i++)
				spin[i] = random.nextDouble() < 0.5 ? 1 : -1;
		}
		
		private int neighborSum(int i) {			
			int x = i % L1;
			int y = i / L1;
			int xs[] = {(x+1)%L1, (x-1+L1)%L1, x, x};
			int ys[] = {y, y, (y+1)%L2, (y-1+L2)%L2};
			
			int acc = 0;
			for (int j = 0; j < 4; j++)
				acc += spin[ys[j]*L1+xs[j]];
			
			return acc;
		}
		
		private void singleStep() {
			int i = random.nextInt(N);
			double dE = 2*spin[i]*neighborSum(i);
			
			if (dE <= 0 || random.nextDouble() < Math.exp(-dE/T))
				spin[i] = -spin[i];
		}
		
		public void step(double mcs) {
			int n = (int) (mcs * N);
			for (int i = 0; i < n; i++)
				singleStep();
			time += mcs;
		}
		
		public NewmanZiff findClusters(int dir) {
			NewmanZiff nz = new NewmanZiff(N);
			for (int i = 0; i < N; i++) {
				if (spin[i] == dir) {
					int x = i % L1;
					int y = i / L1;
					int ip = ((y+1)%L2)*L1 + x;
					if (spin[ip] == dir)
						nz.addBond(i, ip, 0, 1);
					ip = y*L1 + (x+1)%L1;
					if (spin[ip] == dir)
						nz.addBond(i, ip, 1, 0);
				}
			}
			return nz;
		}
	}
}

package chris.Ising;

import chris.util.PrintUtil;
import scikit.graphics.dim2.Grid;
import scikit.jobs.Control;
import scikit.jobs.Job;
import scikit.jobs.Simulation;
import scikit.jobs.params.DoubleValue;

public class isingWolffapp extends Simulation{

	public static final double Tc = 2.0/Math.log(1.0+Math.sqrt(2.0));
								// ~  2.26918531
	public static final int dlength = 10000;
	//public static final int dlength = 100;
    Grid grid = new Grid("Ising spins");
	public int L, N, stack[], spins[];
	public double T, sbar[], data[];
	public boolean draw, cluster[];
	
	
	public static void main(String[] args) {
		new Control(new isingWolffapp(), "Ising Model using Wolff Dynamics");
	}
	
	public void animate() {
		
		T = params.fget("T");
		
		if(draw)
			grid.registerData(L,L,spins);
		
		return;
	}

	public void clear() {
		
		grid.clear();
		return;
	}

	public void load(Control c) {
		params.add("L", 512);
		params.add("T_c");
		params.set("T_c", 2.26918531);
		params.addm("T", new DoubleValue(Tc, 0, 1000));
		params.add("time");
		c.frame(grid);
	}

	public void run() {
		
		L = params.iget("L");
		T = params.fget("T");

		N       = L*L;
		stack   = new int[N];
		spins   = new int[N];
		cluster = new boolean[N];
		sbar    = new double[N];
		data    = new double[dlength];
		
		for (int jj = 0 ; jj < N ; jj++)
			spins[jj] = Math.random() > 0.5 ? 1 : -1;			

		draw = false;
		for (int jj = 0 ; jj < 1000 ; jj++){
		//for (int jj = 0 ; jj < 100 ; jj++){
			step();
			if(jj % 100 == 0) {
				params.set("time",jj-1000);
				Job.animate();
			}
		}
			step();
		

		for (int jj = 0 ; jj < dlength ; jj++){
			step();
			takeData(jj);
			if(jj % 100 == 0){
				params.set("time",jj);
				Job.animate();
			}
		}

		PrintUtil.printVectorToFile("/Users/cserino/Desktop/wolff.txt", data);
		Job.signalStop();
		Job.animate();
		
	}
	
	public void step() {
		// copied from kip.ising.dim2.apps.Wolff2DApp
		int nflipped = 0;
		while (nflipped < N/2) {
			double p = 1 - Math.exp(-2/T);

			stack[0] = (int)(Math.random()*N);
			cluster[stack[0]] = true;
			int stack_len = 1;

			while (stack_len-- > 0) {
				int i = stack[stack_len];
				for (int k = 0; k < 4; k++) {
					int ip = neighbor(i, k);
					if (!cluster[ip] && spins[ip] == spins[i] && Math.random() <= p) {
						cluster[ip] = true;
						stack[stack_len++] = ip;
					}
				}
				spins[i] *= -1;
				cluster[i] = false; // cannot be readded to the cluster since spin is now misaligned
				nflipped++;
			}
		}
		return;
	}
	
	public int neighbor(int i, int k) {
		int y = i/L;
		int x = i%L;
		int yp = (y + (k-1)%2 + L) % L;    // (k-1)%2 == {-1,  0, 1, 0}
		int xp = (x + (k-2)%2 + L) % L;    // (k-2)%2 == { 0, -1, 0, 1}
		return yp*L+xp;
	}
	
	public void takeData(int t){
		double tmp   =  0;
		double omega = 0;
		for (int jj = 0 ; jj < N ; jj++){
			sbar[jj] += spins[jj];
			tmp      += sbar[jj]; 
		}
		tmp /= N;
		for (int jj = 0 ; jj < N ; jj++){
			omega += (sbar[jj]-tmp)*(sbar[jj]-tmp);
		}
		double NN = N;
		double tt = t;
		omega = (NN*tt*tt)/omega;
		data[t] = omega;
		return;
	}
	
	

}

package ranjit.umbrella;

import java.util.Random;

import ranjit.ising.spinblock.SpinBlocks2D;
import scikit.jobs.Control;
import scikit.jobs.Simulation;

public class Ising2D extends Simulation {
	SpinBlocks2D spins;
	int L,R;
	double T;
	Random r=new Random();
	
	public Ising2D() {
		spins = new SpinBlocks2D(L, R);
	}

	public void load(Control c){		
		params.add("T",1.0);
		params.add("L",32);
		params.add("R",8);
	}
	@Override
	public void animate() {
	

	}

	@Override
	public void clear() {
	

	}

	@Override
	public void run() {
		while(true){
			//select spin to flip
//			int x=r.nextInt(L*L);
			//calculate energy cost of spin flip
//			double dE=isingDE(x)+umbrellaDE(x);
			//choose random number
		}
	}

	public double isingDE(int x){
		double dE=0;	
		return dE;
	}
	
	public double umbrellaDE(int x){
		double dE=0;
		return dE;		
	}
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		new Control(new Ising2D(),"Ising 2D Umbrella Sampling");

	}

}

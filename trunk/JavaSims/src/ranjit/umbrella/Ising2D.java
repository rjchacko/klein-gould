package ranjit.umbrella;

import java.util.Random;

import ranjit.ising.spinblock.SpinBlocks2D;
import scikit.jobs.Control;
import scikit.jobs.Simulation;

public class Ising2D extends Simulation {
	SpinBlocks2D spins;
	int L,R;
	double T,h,J;
	double E;
	Random r=new Random();
	int windowMin, windowMax;
	
	public Ising2D() {
		spins = new SpinBlocks2D(L, R);
	}

	public void load(Control c){		
		params.add("T",1.0);
		params.add("h",-0.7);
		params.add("L",32);
		params.add("R",8);
		params.add("window min",992);
		params.add("window max",1024);
	}
	@Override
	public void animate() {
	

	}

	@Override
	public void clear() {
	

	}

	@Override
	public void run() {
		int spinsInRange=(2*R+1)*(2*R+1) - 1;
		J = 4.0 / spinsInRange;
		initializeWindow();
		while(true){
			int x=r.nextInt(L*L);
			double dE=isingDE(x)+umbrellaDE(x);
			if(dE<=0 || Math.exp(-dE/T)<r.nextDouble()) spins.flip(x%L, x/L);
		}
	}
	
	public void initializeWindow(){
		int x, y, j=0;
		double dE;
		do{
			x=j%L;
			y=j/L;
			spins.flip(x,y);
			if(R>1) dE = 2*(h + J*(spins.sumInRange(x,y)-1));
			else dE = 2*(h + J*(spins.get((x-1+L)%L,y)+spins.get((x+1+L)%L,y)+spins.get(x,(y-1+L)%L)+spins.get(x,(y+1+L)%L)));
			E+=dE;
			j++;
		}while(!(spins.netSum>windowMin && spins.netSum<windowMax));
	}

	public double isingDE(int x){
		double dE=0;	
		int spin=spins.get(x%L, x/L);
		dE=2*spin*(h + J*(spins.sumInRange(x%L,x/L)-spin));
		return dE;
	}
	
	public double umbrellaDE(int x){
		double dE=0;
		int dM=-2*spins.get(x%L, x/L);
		if(spins.netSum+dM<windowMin || spins.netSum+dM>windowMax) dE=Double.POSITIVE_INFINITY;
		return dE;		
	}
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		new Control(new Ising2D(),"Ising 2D Umbrella Sampling");

	}

}

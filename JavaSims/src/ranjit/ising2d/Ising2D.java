package ranjit.ising2d;

import java.awt.Color;
import java.util.Random;

import ranjit.ising.spinblock.SpinBlocks2D;
import scikit.dataset.Histogram;
import scikit.graphics.ColorGradient;
import scikit.graphics.ColorPalette;
import scikit.graphics.dim2.Grid;
import scikit.graphics.dim2.Plot;
import scikit.jobs.Control;
import scikit.jobs.Job;
import scikit.jobs.Movies;
import scikit.jobs.Simulation;

public class Ising2D extends Simulation {
	SpinBlocks2D spins;
	Grid grid = new Grid ("Ising Lattice");
	Histogram mag=new Histogram(0.00001);
	Plot magPlot = new Plot("Magnetizations");
	
	int L,R;
	double T,h,J;
	double E;
	Random r=new Random();

	public Ising2D() {

	}

	public void load(Control c){		
		params.add("T",1.0);
		params.add("h",-0.4);
		params.add("L",256);
		params.add("R",1);
		params.add("MCS",10000);
		params.add("mcs");
		c.frame(grid,magPlot);
	}

	public void animate() {
		ColorPalette palette = new ColorPalette();
		palette.setColor(0,Color.BLACK);
		palette.setColor(1,Color.WHITE);

		ColorGradient smooth = new ColorGradient();
		int allSpins[]=spins.blocksAtScale(0);
		for (int i=0 ; i<spins.N ; i++){
			smooth.getColor(allSpins[i],-1,1);
		}
		grid.setColors(smooth);
		grid.registerData(L,L,allSpins);
		magPlot.registerPoints("Magnetization", mag, Color.RED);
	}


	public void clear() {
		grid.clear();
		mag.clear();
		magPlot.clear();
	}


	public void run() {
		R=params.iget("R");
		L=params.iget("L");
		T=params.fget("T");
		h=params.fget("h");
		int totalMCS=params.iget("MCS");

		int spinsInRange=(2*R+1)*(2*R+1) - 1;
		if(R>1)J = 4.0 / spinsInRange;
		else J=1;

		spins = new SpinBlocks2D(L, R);

		for(int mcs=0;mcs<totalMCS;mcs++){
			for(int k=0;k<L*L;k++){
				int nextSpin=r.nextInt(L*L);
				int x=nextSpin%L;
				int y=nextSpin/L;
				
				double dE=isingDE(nextSpin);	
				boolean decreaseEnergy=dE<=0;
				boolean acceptThermal=r.nextDouble()<Math.exp(-dE/T);
				if( decreaseEnergy || acceptThermal)spins.flip(x, y);
				double m=(double)spins.netSum/(L*L);
				mag.accum(m);
				Job.animate();
			}
			params.set("mcs", mcs);			
		}
	}

	public double isingDE(int location){
		double dE=0;	
		int x=location%L;
		int y=location/L;
		int spin=spins.get(x, y);
		if(R>1){
			dE=2*spin*(h + J*(spins.sumInRange(x,y)-spin));
		}
		else {
			dE = 2*spin*(h + J*(spins.get((x-1+L)%L,y)+spins.get((x+1+L)%L,y)+spins.get(x,(y-1+L)%L)+spins.get(x,(y+1+L)%L)));
		}
		return dE;
	}

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		Control c=new Control(new Ising2D(),"Ising 2D Umbrella Sampling");
		new Movies(c);
	}

}

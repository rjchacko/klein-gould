package ranjit.umbrella;

import java.awt.Color;
import java.util.Random;

import ranjit.ising.spinblock.SpinBlocks2D;
import scikit.dataset.Accumulator;
import scikit.dataset.Histogram;
import scikit.graphics.ColorGradient;
import scikit.graphics.ColorPalette;
import scikit.graphics.dim2.Grid;
import scikit.graphics.dim2.Plot;
import scikit.jobs.Control;
import scikit.jobs.Job;
import scikit.jobs.Simulation;

public class Ising2D extends Simulation {
	SpinBlocks2D spins;
	Grid grid = new Grid ("Ising Lattice");
	Histogram mag=new Histogram(1);
    Plot magPlot = new Plot("Magnetizations");
    Accumulator freeEnergy=new Accumulator(1);
    Plot freeEnergyPlot=new Plot("Free Energy");
    
	int L,R;
	double T,h,J;
	double E;
	Random r=new Random();
	int windowMin, windowMax;
	int mcs;
	
	public Ising2D() {
		
	}

	public void load(Control c){		
		params.add("T",1.0);
		params.add("h",-0.7);
		params.add("L",32);
		params.add("R",8);
		params.add("window min",992);
		params.add("window max",1024);
		c.frame(grid,magPlot,freeEnergyPlot);
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
		magPlot.registerPoints("Magnetizations", mag, Color.RED);
		freeEnergyPlot.registerPoints("Free Energy",freeEnergy,Color.RED);
	}

	
	public void clear() {
		grid.clear();
		mag.clear();
		magPlot.clear();
		freeEnergy.clear();
		freeEnergyPlot.clear();
	}

	
	public void run() {
		R=params.iget("R");
		L=params.iget("L");
		T=params.fget("T");
		h=params.fget("h");
		windowMin=params.iget("window min");
		windowMax=params.iget("window max");
		spins = new SpinBlocks2D(L, R);	
		int spinsInRange=(2*R+1)*(2*R+1) - 1;
		J = 4.0 / spinsInRange;
		initializeWindow();
		Job.animate();
		int steps=0;
		while(true){
			int x=r.nextInt(L*L);
			double dE=isingDE(x)+umbrellaDE(x);
			if(dE<=0 || r.nextDouble()<Math.exp(-dE/T)) spins.flip(x%L, x/L);
			mag.accum(spins.netSum);
						
			Job.animate();
			steps++;
			if(steps%(L*L)==0){
				mcs++;
				double data[]=mag.copyData();
				for(int i=0;i<data.length;i+=2){
					freeEnergy.accum(data[i], -Math.log(data[i+1]));
				}
			}
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
		}while(spins.netSum!=windowMax-(windowMax-windowMin)/2.);
		
		for(int cnt=0; cnt<1000; cnt++){
			for(int i=0;i<L*L;i++){
				dE=isingDE(x)+umbrellaDE(x);
				if(dE<=0 || r.nextDouble()<Math.exp(-dE/T)) spins.flip(x%L, x/L);		
			}
		}
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

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
	Histogram mag[];
    Plot magPlot = new Plot("Magnetizations");
    Accumulator freeEnergy[];
    Plot freeEnergyPlot=new Plot("Free Energy");
    
	int L,R;
	double T,h,J;
	double E;
	Random r=new Random();
	int windowMin, windowMax, currentWindow;
	
	public Ising2D() {
		
	}

	public void load(Control c){		
		params.add("T",1.0);
		params.add("h",-0.7);
		params.add("L",32);
		params.add("R",1);
		params.add("number of windows", 256);
		params.add("window width",8);
		params.add("MCS per window", 1000);
		params.add("mcs");
		params.add("current window");
		params.add("prefix","/Users/rjchacko/Desktop/data/");
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
		if(currentWindow%2==0)magPlot.registerPoints("Magnetizations"+ currentWindow, mag[currentWindow], Color.RED);
		else magPlot.registerPoints("Magnetizations"+ currentWindow, mag[currentWindow], Color.BLUE);
		if(currentWindow%2==0)freeEnergyPlot.registerPoints("Free Energy"+currentWindow,freeEnergy[currentWindow],Color.RED);
		else freeEnergyPlot.registerPoints("Free Energy"+currentWindow,freeEnergy[currentWindow],Color.BLUE);
	}

	
	public void clear() {
		grid.clear();
		magPlot.clear();
		freeEnergyPlot.clear();
	}

	
	public void run() {
		R=params.iget("R");
		L=params.iget("L");
		T=params.fget("T");
		h=params.fget("h");
		int numberOfWindows=params.iget("number of windows");
		int mcsPerWindow=params.iget("MCS per window");
		String prefix=params.sget("prefix");
		int windowSpacing=2*L*L/numberOfWindows;
		int windowWidth=params.iget("window width");
		
		int spinsInRange=(2*R+1)*(2*R+1) - 1;
		if(R>1)J = 4.0 / spinsInRange;
		else J=1;
		freeEnergy=new Accumulator[numberOfWindows];
		for(int l=0;l<numberOfWindows;l++){
			freeEnergy[l]=new Accumulator(1);
		}
		mag=new Histogram[numberOfWindows];
		for(int l=0;l<numberOfWindows;l++){
			mag[l]=new Histogram(1);
		}
		for(int j=0;j<numberOfWindows;j++){
			params.set("current window", j);
			initializeWindow(windowSpacing,windowWidth, j);
			Job.animate();
			for(int mcs=0;mcs<mcsPerWindow;mcs++){
				for(int k=0;k<L*L;k++){
					int nextSpin=r.nextInt(L*L);
					double dE=isingDE(nextSpin)+umbrellaDE(nextSpin);	
					boolean decreaseEnergy=dE<=0;
					boolean acceptThermal=r.nextDouble()<Math.exp(-dE/T);
					int x=nextSpin%L;
					int y=nextSpin/L;
					if( decreaseEnergy || acceptThermal) spins.flip(x, y);
					mag[currentWindow].accum(spins.netSum);		
					Job.animate();
				}
				params.set("mcs", mcs);	
				double data[]=mag[currentWindow].copyData();
				for(int i=0;i<data.length;i+=2){
					freeEnergy[currentWindow].accum(data[i], -T*Math.log(data[i+1]));
				}
			}
			
			freeEnergyPlot.saveDataset2(freeEnergy[currentWindow], prefix+"Free Energy"+currentWindow+".txt");
			magPlot.saveDataset2(mag[currentWindow], prefix+"Magnetization"+currentWindow+".txt");
		}
		
	}
	
	public void initializeWindow(int windowSpacing, int windowWidth, int windowNumber){
		int x, y, j=0;
		double dE;
		currentWindow=windowNumber;
		int phi0=L*L-(windowNumber+1)*windowSpacing;
		windowMax=phi0+windowWidth;
		windowMin=phi0-windowWidth;
		spins = new SpinBlocks2D(L, R);
		do{
			x=j%L;
			y=j/L;
			spins.flip(x,y);
			if(R>1) dE = 2*(h + J*(spins.sumInRange(x,y)-1));
			else dE = 2*(h + J*(spins.get((x-1+L)%L,y)+spins.get((x+1+L)%L,y)+spins.get(x,(y-1+L)%L)+spins.get(x,(y+1+L)%L)));
			E+=dE;
			j++;
		}while(spins.netSum!=windowMax-(windowMax-windowMin)/2.);
		
		for(int cnt=0; cnt<500; cnt++){
			for(int i=0;i<L*L;i++){
				int nextSpin=r.nextInt(L);
				dE=isingDE(nextSpin)+umbrellaDE(nextSpin);
				boolean decreaseEnergy=dE<=0;
				boolean acceptThermal=r.nextDouble()<Math.exp(-dE/T);
				x=nextSpin%L;
				y=nextSpin/L;
				if( decreaseEnergy || acceptThermal) spins.flip(x, y);		
			}
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
			dE = 2*spin*(h + J*(spins.get(x,y)+spins.get((x+1+L)%L,y)+spins.get(x,(y-1+L)%L)+spins.get(x,(y+1+L)%L)));
		}
		return dE;
	}
	
	public double umbrellaDE(int location){
		double dE=0;
		int x=location%L;
		int y=location/L;
		int dM=-2*spins.get(x, y);
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

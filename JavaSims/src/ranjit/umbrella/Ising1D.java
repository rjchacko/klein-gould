package ranjit.umbrella;

import java.awt.Color;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Random;

import ranjit.ising.spinblock.SpinBlocks1D;
import scikit.dataset.Accumulator;
import scikit.dataset.DataSet;
import scikit.dataset.DatasetBuffer;
import scikit.dataset.Histogram;
import scikit.graphics.dim2.Plot;
import scikit.jobs.Control;
import scikit.jobs.Job;
import scikit.jobs.Simulation;
import scikit.util.FileUtil;

public class Ising1D extends Simulation {
	SpinBlocks1D spins;
	Histogram mag[];
	Plot magPlot = new Plot("Magnetizations");
	Accumulator freeEnergy[];
	Plot freeEnergyPlot=new Plot("Free Energy");
	boolean randomic=true;

	int L,R;
	double T,h,J;
	double E;
	Random r=new Random();
	int windowMin, windowMax, currentWindow;
	int firstWindow, lastWindow;
	public Ising1D() {

	}

	public void load(Control c){		
		params.add("T",1.0);
		params.add("h",-0.4);
		params.add("L",256);
		params.add("R",1);
		params.add("number of windows", 32768);
		params.add("first window",0);
		params.add("last window", 100);
		params.add("window width",4);
		params.add("MCS per window", 1000);
		params.add("mcs");
		params.add("current window");
		params.add("phi0");
		params.add("prefix","/Users/rjchacko/Desktop/data/");
		c.frame(magPlot,freeEnergyPlot);
	}

	public void animate() {
		if(currentWindow%2==0)magPlot.registerPoints("Magnetizations"+ currentWindow, mag[currentWindow], Color.RED);
		else magPlot.registerPoints("Magnetizations"+ currentWindow, mag[currentWindow], Color.BLUE);
		if(currentWindow%2==0)freeEnergyPlot.registerPoints("Free Energy"+currentWindow,freeEnergy[currentWindow],Color.RED);
		else freeEnergyPlot.registerPoints("Free Energy"+currentWindow,freeEnergy[currentWindow],Color.BLUE);
	}


	public void clear() {
		magPlot.clear();
		freeEnergyPlot.clear();
	}

	private void saveDataset(DataSet data, String fname) {
		try {
			if (fname != null) {
				PrintWriter pw = FileUtil.pwFromString(fname);
				FileUtil.writeColumns(pw, data.copyData().columns());
				pw.close();
			}
		} catch (IOException e) {}
	}


	public void run() {
		R=params.iget("R");
		L=params.iget("L");
		T=params.fget("T");
		h=params.fget("h");
		int numberOfWindows=params.iget("number of windows");
		int mcsPerWindow=params.iget("MCS per window");
		String prefix=params.sget("prefix");
		int windowSpacing=2*L/numberOfWindows;
		int windowWidth=params.iget("window width");
		int firstWindow=params.iget("first window");
		int lastWindow=params.iget("last window");

		int spinsInRange=2*R;
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
		int init[]=new int[L];
		spins = new SpinBlocks1D(L, R, 1);

		do{
			int x=r.nextInt(L);
			if(spins.get(x)==1)spins.flip(x);
		}while(spins.netSum>L-windowSpacing);
		
		init=spins.getAll().clone();

		for(int j=firstWindow;j<=lastWindow;j++){
			params.set("current window", j);
			initializeWindow(windowSpacing,windowWidth, j, init);

			int phi0=L-(j+1)*windowSpacing;
			params.set("phi0", phi0);
			Job.animate();
			boolean saved=false;
			for(int mcs=0;mcs<mcsPerWindow;mcs++){
				for(int k=0;k<L;k++){
					int nextSpin=r.nextInt(L);
					double dE=isingDE(nextSpin)+umbrellaDE(nextSpin);	
					boolean decreaseEnergy=dE<=0;
					boolean acceptThermal=r.nextDouble()<Math.exp(-dE/T);
					if( decreaseEnergy || acceptThermal)spins.flip(nextSpin);
					mag[currentWindow].accum(spins.netSum);

					if(!randomic && spins.netSum<phi0 && !saved && mcs>100){
						init=spins.getAll().clone();
						saved=true;
					}
					Job.animate();
				}
				params.set("mcs", mcs);	
				DatasetBuffer data=mag[currentWindow].copyData();
				for(int i=0;i<data.size();i++){
					freeEnergy[currentWindow].accum(data.x(i), -T*Math.log(data.y(i)));
				}
			}

			saveDataset(freeEnergy[currentWindow], prefix+"Free Energy"+currentWindow+".txt");
			saveDataset(mag[currentWindow], prefix+"Magnetization"+currentWindow+".txt");
		}

	}

	public void initializeWindow(int windowSpacing, int windowWidth, int windowNumber,int init[]){
		double dE;
		currentWindow=windowNumber;
		int phi0=L-(windowNumber+1)*windowSpacing;
		windowMax=phi0+windowWidth;
		windowMin=phi0-windowWidth;
		spins = new SpinBlocks1D(L, R, 1);

		if(!randomic){
			for(int i=0;i<init.length;i++){
				if(init[i]==-1)spins.flip(i);
			}
		}
		else{
			while(spins.netSum>phi0){
				int spin=r.nextInt(L);
				if(spins.get(spin)==1) spins.flip(spin);
			}
		}

		for(int cnt=0; cnt<100; cnt++){
			for(int i=0;i<L;i++){
				int nextSpin=r.nextInt(L);
				dE=isingDE(nextSpin)+umbrellaDE(nextSpin);
				boolean decreaseEnergy=dE<=0;
				boolean acceptThermal=r.nextDouble()<Math.exp(-dE/T);
				if( decreaseEnergy || acceptThermal) spins.flip(nextSpin);		
			}
		}
	}

	public double isingDE(int location){
		double dE=0;	
		int spin=spins.get(location);
		
		dE=2*spin*(h + J*(spins.sumInRange(location)-spin));
		return dE;
	}

	public double umbrellaDE(int location){
		double dE=0;
		int dM=-2*spins.get(location);
		if(spins.netSum+dM<windowMin || spins.netSum+dM>windowMax) dE=Double.POSITIVE_INFINITY;
		return dE;		
	}
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		new Control(new Ising1D(),"Ising 1D Umbrella Sampling");

	}

}

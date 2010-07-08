package kang.umbrella;

import java.awt.Color;
import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Random;


import javax.imageio.ImageIO;

import kang.umbrella.spinblock.SpinBlocks2D;
import scikit.dataset.Accumulator;
import scikit.dataset.DataSet;
import scikit.dataset.DatasetBuffer;
import scikit.dataset.Histogram;
import scikit.graphics.ColorGradient;
import scikit.graphics.ColorPalette;
import scikit.graphics.dim2.Grid;
import scikit.graphics.dim2.Plot;
import scikit.jobs.Control;
import scikit.jobs.Job;
import scikit.jobs.Movies;
import scikit.jobs.Simulation;
import scikit.util.FileUtil;

public class Ising2D extends Simulation {
	SpinBlocks2D spins;
	Grid grid = new Grid ("Ising Lattice");
	Histogram mag[];
	Plot magPlot = new Plot("Magnetizations");
	Accumulator freeEnergy[];
	Plot freeEnergyPlot=new Plot("Free Energy");
	boolean randomic=true;

	int L,R;
	double T,h,J;
	double E;
	double dilution; //dilution
	Random r=new Random();
	int windowMin, windowMax, currentWindow;
	int firstWindow, lastWindow;
	public Ising2D() {

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
		params.add("prefix","/Users/liukang2002507/Desktop/umbrelladata/");
		c.frame(grid,magPlot,freeEnergyPlot);
		params.add("dilution",0.1);
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
		dilution=params.fget("dilution");
		
		int numberOfWindows=params.iget("number of windows");
		int mcsPerWindow=params.iget("MCS per window");
		String prefix=params.sget("prefix");
		int windowSpacing=2*L*L/numberOfWindows;
		int windowWidth=params.iget("window width");
		int firstWindow=params.iget("first window");
		int lastWindow=params.iget("last window");

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
		int init[]=new int[L*L];
		spins = new SpinBlocks2D(L, R, dilution,1);

		do{
			int z=r.nextInt(L*L);
			int x=z%L;
			int y=z/L;
			if(spins.get(x, y)==1)spins.flip(x, y);
		}while(spins.netSum>L*L-windowSpacing);
		init=spins.blocksAtScale(0).clone();

		for(int j=firstWindow;j<=lastWindow;j++){
			params.set("current window", j);
			initializeWindow(windowSpacing,windowWidth, j, init);

			int phi0=L*L-(j+1)*windowSpacing;
			params.set("phi0", phi0);
			Job.animate();
			boolean saved=false;
			for(int mcs=0;mcs<mcsPerWindow;mcs++){
				for(int k=0;k<L*L;k++){
					int nextSpin=r.nextInt(L*L);
					double dE=isingDE(nextSpin)+umbrellaDE(nextSpin);	
					boolean decreaseEnergy=dE<=0;
					boolean acceptThermal=r.nextDouble()<Math.exp(-dE/T);
					int x=nextSpin%L;
					int y=nextSpin/L;
					if( decreaseEnergy || acceptThermal)spins.flip(x, y);
					mag[currentWindow].accum(spins.netSum);

					if(!randomic && spins.netSum<phi0 && !saved && mcs>100){
						init=spins.blocksAtScale(0).clone();
						saved=true;
					}
					
					Job.animate();
				}
				params.set("mcs", mcs);			
			}
			DatasetBuffer data=mag[currentWindow].copyData();
			for(int i=0;i<data.size();i++){
				freeEnergy[currentWindow].accum(data.x(i), -T*Math.log(data.y(i)));
			}
			saveDataset(freeEnergy[currentWindow], prefix+"Free Energy"+currentWindow+".txt");
			saveDataset(mag[currentWindow], prefix+"Magnetization"+currentWindow+".txt");
				
			try {
					ImageIO.write(grid.getImage(), "png", new File(prefix+currentWindow+".png"));
			}catch (IOException e) {}
		}

	}

	public void initializeWindow(int windowSpacing, int windowWidth, int windowNumber,int init[]){
		int x, y;
		double dE;
		currentWindow=windowNumber;
		int phi0=L*L-(windowNumber+1)*windowSpacing;
		windowMax=phi0+windowWidth;
		windowMin=phi0-windowWidth;
		spins = new SpinBlocks2D(L, R, dilution,1);

		if(!randomic){
			for(int i=0;i<init.length;i++){
				x=i%L;
				y=i/L;
				if(init[i]==-1)spins.flip(x, y);
			}
		}
		else{
			while(spins.netSum>phi0){
				int spin=r.nextInt(L*L);
				x=spin%L;
				y=spin/L;
				if(spins.get(x, y)==1) spins.flip(x,y);
			}
		}

		for(int cnt=0; cnt<100; cnt++){
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
			dE = 2*spin*(h + J*(spins.get((x-1+L)%L,y)+spins.get((x+1+L)%L,y)+spins.get(x,(y-1+L)%L)+spins.get(x,(y+1+L)%L)));
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
		Control c=new Control(new Ising2D(),"Ising 2D Umbrella Sampling");
		new Movies(c);
	}

}

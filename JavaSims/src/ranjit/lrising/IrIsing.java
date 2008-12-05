package ranjit.lrising;

import java.awt.Color;
import java.io.FileOutputStream;
import java.io.ObjectOutputStream;
import java.util.ArrayList;
import java.util.Random;

import scikit.dataset.Accumulator;
import scikit.dataset.Histogram;
import scikit.graphics.ColorGradient;
import scikit.graphics.ColorPalette;
import scikit.graphics.dim2.Grid;
import scikit.graphics.dim2.Plot;
import scikit.jobs.Control;
import scikit.jobs.Job;
import scikit.jobs.Simulation;

public class IrIsing extends Simulation {

	int spins[]=null;
	int L;
	int N;
	double T;
	double h;
	int M;
	double E;
	double J;
	Grid grid = new Grid ("Ising Lattice");
	Histogram clusterHist=new Histogram(1);
	
	Accumulator clustTotalAccum= new Accumulator();
	Accumulator magAccum=new Accumulator();
	Accumulator momentAccum=new Accumulator();
	Histogram clustBins[]=new Histogram[60];
	Accumulator largestClustAccum=new Accumulator();
	
	Plot magPlot = new Plot("Magnetizations");
	Plot momentPlot=new Plot("Second Moment");
	Plot clustTotalPlot=new Plot("Total number of clusters");
	Plot largestClustPlot=new Plot("Mass of largest cluster");
	Plot clustPlot=new Plot("Cluster Numbers");
	
	lrClusters clusters;
	ArrayList<int[]> history=new ArrayList<int[]>();
	
	public void animate() {
		ColorPalette palette = new ColorPalette();
		palette.setColor(0,Color.BLACK);
		palette.setColor(1,Color.WHITE);

		ColorGradient smooth = new ColorGradient();
		for (int i=0 ; i<spins.length; i++){
			smooth.getColor(spins[i],-1,1);
		}
		grid.setColors(smooth);
		grid.registerData(L,L,spins);
		magPlot.registerPoints("magnetization", magAccum, Color.red);
		momentPlot.registerPoints("moment", momentAccum, Color.blue);
		clustTotalPlot.registerPoints("allclusters", clustTotalAccum, Color.green);
		largestClustPlot.registerPoints("largestclust", largestClustAccum,Color.red);
		for(int i=0;i<60;i++){
			clustPlot.registerPoints("numClust"+i, clustBins[i], new Color((int)Math.exp(2*i)));
		}
	}

	@Override
	public void clear() {
		grid.clear();
		magAccum.clear();
		magPlot.clear();
		clusterHist.clear();
		momentAccum.clear();
		momentPlot.clear();
		clustTotalPlot.clear();
		clustTotalAccum.clear();
		
		for(int i=0;i<L*L;i++){
			spins[i]=1;
		}
		
	}

	@Override
	public void load(Control c) {
		params.add("T",1.78);
		params.add("h",-1.2);
		params.add("L",32);
		params.add("MCS",100000);
		params.add("mcs");
		params.add("M");
		params.add("Bond probability");
		
		c.frameTogether("graphs",grid,magPlot,momentPlot,clustTotalPlot,largestClustPlot,clustPlot);
	}

	@Override
	public void run() {
		L=params.iget("L");
		spins=new int[L*L];
		J=4./(L*L);
		T=params.fget("T");
		h=params.fget("h");
		int totalMCS=params.iget("MCS");
		Random r=new Random();
		M=L*L;
		E=-0.5*J*M*M-h*M;
		int dM=-2;
		for(int i=0;i<L*L;i++){
			spins[i]=1;
		}
		double m=(double)M/(double)(L*L);
		double x=1.+m;
		double bondProbability=1.0-Math.exp(-2.*(J/T)*x);
		
		for(int i=0;i<60;i++){
			clustBins[i]=new Histogram(1);
		}
			
		clusters=new lrClusters(L,bondProbability);
		for(int mcs=0;mcs<totalMCS;mcs++){
			
			for(int k=0;k<L*L;k++){
				int nextSpin=r.nextInt(L*L);
				int sign=spins[nextSpin];
				int newM=M+sign*dM;
				double newEnergy=(-0.5*J*newM*newM-h*newM);
				double dE=newEnergy-E;
				boolean decreaseEnergy=dE<=0;
				boolean acceptThermal=r.nextDouble()<Math.exp(-dE/T);
				if( decreaseEnergy || acceptThermal){
					spins[nextSpin]*=-1;
					M=newM;
					E=newEnergy;
				}
			}
			params.set("mcs", mcs);
			params.set("M", M);
			history.add(spins.clone());
			if(history.size()==100)history.remove(0);
			//save spin configurations for intervention
			try {
	            FileOutputStream fos = new FileOutputStream("history");
	            ObjectOutputStream oos = new ObjectOutputStream(fos);
	            oos.writeObject(history);
	            oos.flush();
	            fos.close();
			}
			catch (Throwable e) {
	            System.err.println("exception thrown");
			}	
			
			m=(double)M/(double)(L*L);
			x=1.+m;
			bondProbability=1.0-Math.exp(-2.*(J/T)*x);
			clusters.bondProbability=bondProbability;
			params.set("Bond probability", bondProbability);
		
			clusters.newLattice();
			clusterHist.clear();
			for(int i=0;i<spins.length;i++){
				if(spins[i]==-1) {
					clusters.addSite(i);
				}
			}
			
			int allclusters=0;
			
			for(int i=1;i<clusters.numClusters.length;i++){
				clustBins[i%10].accum(mcs, clusters.numClusters[i]);
				allclusters+=clusters.numClusters[i];
			}
			magAccum.accum(mcs, m);
			momentAccum.accum(mcs, Math.sqrt(clusters.secondClusterMoment));
			clustTotalAccum.accum(mcs, allclusters);
			largestClustAccum.accum(mcs, clusters.massLargestCluster);
			Job.animate();
		}

	}


	public static void main(String[] args) {
		new Control(new IrIsing(),"Ising 2D Umbrella Sampling");
	}

}

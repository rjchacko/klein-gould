package ranjit.cwising;

import java.awt.Color;
import java.util.ArrayList;

import scikit.dataset.Histogram;

import scikit.graphics.dim2.Plot;
import scikit.jobs.*;

public class cwIsingFC extends Simulation {

	int spins[]=null;
	int L;
	int N;
	double T;
	double h;
	int M;
	double E;
	double J;
	
	Histogram clusterHist=new Histogram(1);
	Plot clustPlot=new Plot("Cluster Size Histogram");
	
	cwClusters clusters;
	ArrayList<int[]> history=new ArrayList<int[]>();
	
	public void animate() {
		clustPlot.registerBars("numclust", clusterHist, Color.red);
	}

	@Override
	public void clear() {
		clusterHist.clear();
		clustPlot.clear();
		for(int i=0;i<L*L;i++){
			spins[i]=1;
		}
		
	}

	@Override
	public void load(Control c) {
		params.add("T",1.78);
		params.add("h",-1.22);
		params.add("L",256);
		params.add("M0", 42475);
		params.add("M");
		params.add("Bond probability");
		c.frame(clustPlot);
	}

	public void run() {
		L=params.iget("L");
		spins=new int[L*L];
		J=4./(L*L);
		T=params.fget("T");
		h=params.fget("h");
		int M0=params.iget("M0");
		for(int i=0;i<spins.length;i++){
			spins[i]=1;
		}
		M=L*L;
		params.set("M", M);
		int j=0;
		while(M>M0){
			spins[j]*=-1;
			M-=2;
			params.set("M", M);
			j++;
		}
		
		double m=(double)M/(double)(L*L);
		double x=1.+m;
		double bondProbability=1.0-Math.exp(-2.*(J/T)*x);
		clusters=new cwClusters(L,bondProbability);
		clusters.bondProbability=bondProbability;
		params.set("Bond probability", bondProbability);
		clusters.newLattice();
		clusterHist.clear();
		for(int i=0;i<spins.length;i++){
			if(spins[i]==-1) {
				clusters.addSite(i);
			}
		}
			
		for(int i=1;i<clusters.numClusters.length;i++){
			clusterHist.accum(i,clusters.numClusters[i]);			
		}
		Job.animate();
		

	}


	public static void main(String[] args) {
		Control c=new Control(new cwIsingFC(),"Curie Weiss Ising");
		new Movies(c);
	}

}

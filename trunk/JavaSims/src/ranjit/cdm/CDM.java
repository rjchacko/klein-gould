package ranjit.cdm;
import java.awt.Color;
import java.util.Random;

import scikit.dataset.Accumulator;
import scikit.dataset.Histogram;
import scikit.graphics.dim2.Plot;
import scikit.jobs.Control;
import scikit.jobs.Job;
import scikit.jobs.Simulation;

public class CDM extends Simulation{
	Histogram magnetization=new Histogram(1);
	Accumulator N=new Accumulator(1);
	Plot numDown=new Plot("Number of spins in stable state");
	Plot mag=new Plot("magnetization");
	double Nl[];
	double h;
	double T;
	int dimension;
	
	public static void main(String[] args) {
		new Control(new CDM(), "Classical Droplet Model");
	}
	public void load(Control c){
		params.add("h",-0.05);
		params.add("T",1.0);
		params.add("dimension",2);
		c.frame(numDown);
		c.frame(mag);
	}
	
	@Override
	public void animate() {
		numDown.registerPoints("Number of Down Spins", N, Color.RED);	
		mag.registerBars("magnetization", magnetization, Color.RED);
	}

	@Override
	public void clear() {
		numDown.clear();
		magnetization.clear();
	}

	@Override
	public void run() {
		Random r= new Random();
		h=params.fget("h");
		T=params.fget("T");
		dimension=params.iget("dimension");
		double maxL=1.25*criticalLength(h);
		System.out.println("maxL="+maxL);
		double dE;
		Nl=new double[(int) maxL];
		numDown.setAutoScale(true);
		int dN;
		int nminus=0;
		while(true){
			
			for(int i=0;i<Nl.length;i++){
				if(r.nextDouble()<0.5) dN=1;
				else dN=-1;
				dE=deltaE(i);
				if(Math.exp(-dE/T)<r.nextDouble()){
					if(Nl[i]>0){
						Nl[i]+=dN;
						nminus+=i*dN;
					}
					else if(Nl[i]==0 && dN>0){
						Nl[i]+=dN;
						nminus+=i*dN;
					}
				}
				N.accum(i, Nl[i]);
			}	
			
			magnetization.accum(100-2*nminus);
			Job.animate();
			N.clear();
		}
	}
	
	double criticalLength(double h){
		int d=dimension;
		return Math.pow((d-1)/(2*d*Math.abs(h)),d);		
	}
	
	double deltaE(int l){		
		return (h*l+Math.pow(l, (dimension-1)/dimension));
	}

}

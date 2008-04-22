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
	Histogram magnetization=new Histogram(10);
	Accumulator N=new Accumulator(1);
	Accumulator mt=new Accumulator(1);
	Accumulator delta=new Accumulator(1);
	Plot numDown=new Plot("Number of spins in stable state");
	Plot mag=new Plot("magnetization");
	Plot moft=new Plot("m(t)");
	Plot delE=new Plot("delta E");
	int Nl[];
	double deltaE[];
	double h;
	double T;
	int dimension;
	int total;
	int mcs=0;
	public static void main(String[] args) {
		new Control(new CDM(), "Classical Droplet Model");
	}
	public void load(Control c){
		params.add("h",-0.1);
		params.add("T",10.0);
		params.add("L",10000);
		params.add("dimension",2);
		params.add("max L");
		c.frameTogether("Classical Droplet Model", delE, moft, mag, numDown);
		numDown.setAutoScale(false);
	}
	
	@Override
	public void animate() {
		params.set("max L", Math.pow(Math.abs(h), -dimension));
		numDown.registerPoints("Number of Down Spins", N, Color.RED);	
		mag.registerBars("magnetization", magnetization, Color.RED);
		moft.registerPoints("m(t)", mt, Color.RED);
		delE.registerPoints("delta E", delta, Color.RED);
	}

	@Override
	public void clear() {
		numDown.clear();
		magnetization.clear();
		delE.clear();
		delta.clear();
		mag.clear();
		moft.clear();
		N.clear();
		mt.clear();
		mcs=0;
	}

	@Override
	public void run() {
		Random r= new Random();
		h=params.fget("h");
		T=params.fget("T");
		dimension=params.iget("dimension");
		total=(int) Math.pow(params.iget("L"), dimension);
		
//		double maxL= criticalLength();
		double maxL= Math.pow(Math.abs(h), -dimension)*2;
		
		
		Nl=new int[(int)maxL];
		
		deltaE=new double[(int)maxL];
		for(int i=0;i<deltaE.length;i++){
			deltaE[i]=deltaE(i);
			delta.accum(i, deltaE[i]);
		}
		Job.animate();
		
		int mag=total;
		mag = equilibrate(r, mag);
		while(true){
			mag = oneMCS(r, mag);			
			magnetization.accum(mag);
			Job.animate();
			N.clear();
		}
	}
	
	private int equilibrate(Random r, int mag) {
		int dN;
		for(int j=0;j<100000;j++){
			for(int i=0;i<Nl.length;i++){
				int index=r.nextInt(Nl.length);
				if(r.nextDouble()<0.5)dN=1; else dN=-1;			
				int dM=-2*index*dN;
				double dE=dN*deltaE[index];
				int newN=Nl[index]+dN;
				int newM=mag+dM;
				if((dE<0 ||(r.nextDouble()<Math.exp(-dE/T))) && newN>0 && newM>total*0.9){
					Nl[index]+=dN;
					mag=newM;
				}
			}
		}
		return mag;
	}
	
	private int oneMCS(Random r, int mag) {
		int dN;
		for(int i=0;i<Nl.length;i++){
			int index=r.nextInt(Nl.length);
			if(r.nextDouble()<0.5)dN=1; else dN=-1;			
			int dM=-2*index*dN;
			double dE=dN*deltaE[index];
			int newN=Nl[index]+dN;
			int newM=mag+dM;
			if((dE<0 ||(r.nextDouble()<Math.exp(-dE/T))) && newN>0 && newM>total*0.9){
				Nl[index]+=dN;
				mag=newM;
			}
			N.accum(index, Nl[index]);
			mt.accum(mcs,mag/(double)total);
			mcs++;
		}
		return mag;
	}
	
	double criticalLength(){
		int d=dimension;
		return Math.pow((d-1)/(2*d*Math.abs(h)),d);		
	}
	
	double deltaE(int l){	
		double x=((double)dimension-1)/(double)dimension;
		return h*l+Math.pow(l, x);
	}

}

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
	Random r= new Random();
	Histogram magnetization=new Histogram(0.001);
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
	double totalL;
	double volume;
	int mcs=0;
	public static void main(String[] args) {
		new Control(new CDM(), "Classical Droplet Model");
	}
	public void load(Control c){
		params.add("h",0.1);
		params.add("T",1.0);
		params.add("L",100);
		params.add("dimension",2);
		c.frameTogether("Classical Droplet Model", delE, moft, mag, numDown);
		numDown.setAutoScale(false);
	}
	
	@Override
	public void run() {
		
		initialize();
		Job.animate();	
		
		while(true){
			oneMCS();		
			
			Job.animate();
			N.clear();
		}
	}
	
	@Override
	public void animate() {
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
	
	private void initialize() {
		h=params.fget("h");
		T=params.fget("T");
		dimension=params.iget("dimension");
		volume=(int) Math.pow(params.iget("L"), dimension);
		totalL=0;
		
		if(h<0){
			double maxL= criticalLength();
			Nl=new int[(int)maxL];
			deltaE=new double[(int)maxL];
		}
		else{
			Nl=new int[(int)volume];
			deltaE=new double[(int)volume];
		}
		
		for(int i=0;i<deltaE.length;i++){
			deltaE[i]=deltaE(i);
			delta.accum(i, deltaE[i]);
		}
	}

	private void oneMCS() {
		int dN;
		for(int i=0;i<Nl.length;i++){
			int l=r.nextInt(Nl.length);
			if(r.nextDouble()<0.5)dN=1; else dN=-1;			
			
			double dE=dN*deltaE[l];
			
			if(Nl[l]+dN>=0 && totalL+dN*l<volume){				
				if((dE<0 ||(r.nextDouble()<Math.exp(-dE/T)))){
					Nl[l]+=dN;
					totalL+=dN*l;
				}
			}
			N.accum(l, Nl[l]);
			
			double m=1-2*totalL/volume;
			
			mt.accum(mcs,m);
			magnetization.accum(m);
			mcs++;
		}
		return;
	}
	
	double criticalLength(){
		double d=dimension;
		return Math.pow((d-1)/(2*d*Math.abs(h)),d);		
	}
	
	double deltaE(int l){	
		double x=((double)dimension-1)/(double)dimension;
		return h*l+Math.pow(l, x);
	}

}

package ranjit.cdm;

import java.awt.Color;
import java.util.Random;

import scikit.dataset.Accumulator;
import scikit.dataset.Histogram;
import scikit.graphics.dim2.Plot;
import scikit.jobs.Control;
import scikit.jobs.Job;
import scikit.jobs.Simulation;

public class cdmXi extends Simulation {
	Random r= new Random();
	Histogram magnetization=new Histogram(0.001);
	Accumulator N=new Accumulator(1);
	Accumulator mt=new Accumulator(1);
	Accumulator delta=new Accumulator(1);
	Accumulator xi=new Accumulator(0.00001);
	
	Plot numDown=new Plot("Number of spins in stable state");
	Plot mag=new Plot("magnetization");
	Plot moft=new Plot("m(t)");
	Plot XI=new Plot("Xi(h)");
	
	int Nl[];
	double deltaE[];
	double h;
	double T;
	int dimension;
	int runLength;
	double totalL;
	double volume;
	double l2accum, laccum;
	int mcs=0;
	
	public static void main(String[] args) {
		new Control(new cdmXi(),"Susceptibility of CDM");
	}
	
	public void load(Control c){
		params.add("T",1.0);
		params.add("L",100);
		params.add("dimension",2);
		params.add("MCS per measurement", 1000);
		params.add("h");
		c.frame(moft, mag,XI);
		numDown.setAutoScale(false);
	}
	
	@Override
	public void run() {
		
		initialize();
		Job.animate();	
		
		for(h=0;h>-0.05;h-=0.001){
			params.set("h", h);
			clear();
			initialize();
			for(int i=0;i<runLength;i++){
				oneMCS();
				Job.animate();
				N.clear();			
			}
			double susceptibility=l2accum/runLength-(laccum/runLength)*(laccum/runLength);
			xi.accum(h,susceptibility);
			
		}
	}
	
	@Override
	public void animate() {	
		numDown.registerPoints("Number of Down Spins", N, Color.RED);	
		mag.registerBars("magnetization", magnetization, Color.RED);
		moft.registerPoints("m(t)", mt, Color.RED);
		XI.registerPoints("Xi(h)",xi, Color.RED);
	}

	@Override
	public void clear() {
		numDown.clear();
		magnetization.clear();
		delta.clear();
		mag.clear();
		moft.clear();
		N.clear();
		mt.clear();
		mcs=0;
	}
	
	private void initialize() {
		l2accum=0;
		laccum=0;
	
		T=params.fget("T");
		dimension=params.iget("dimension");
		runLength=params.iget("MCS per measurement");
		
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
		}
		mcs++;
		double m=1-2*totalL/volume;
		laccum+=totalL;
		l2accum+=totalL*totalL;
		
		mt.accum(mcs,m);
		magnetization.accum(m);
		
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

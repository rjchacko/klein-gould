package ranjit.cdm;

import java.awt.Color;

import scikit.dataset.Accumulator;
import scikit.graphics.dim2.Plot;
import scikit.jobs.Control;
import scikit.jobs.Simulation;
import scikit.jobs.Job;

public class cdmFreeEnergy extends Simulation {
	Accumulator f=new Accumulator(0.00001);
	Plot F=new Plot("f(h)");
	Accumulator fprime=new Accumulator(0.00001);
	Plot Fprime=new Plot("f'(h)");
	Accumulator clength=new Accumulator(0.00001);
	Plot clPlot=new Plot("Critical length");
	Accumulator g=new Accumulator(0.00001);
	Plot G=new Plot("g(l)");
	
	public void load(Control c){
		params.add("T",1.0);
		params.add("hmin", -0.1);
		params.add("hmax", 0.1);
		params.add("volume", 1000);
		params.add("dimension",2);
		c.frame(F,Fprime,clPlot,G);
	}
	
	@Override
	public void animate() {
		F.registerPoints("Free Energy", f, Color.RED);
		Fprime.registerPoints("dF/dH", fprime, Color.RED);
		clPlot.registerPoints("", clength, Color.RED);
		G.registerPoints("g(l)", g, Color.RED);
	}

	@Override
	public void clear() {
		f.clear();
		F.clear();
		fprime.clear();
		Fprime.clear();
		clength.clear();
		clPlot.clear();
		g.clear();
		G.clear();		
	}

	@Override
	public void run() {
		double T=params.fget("T");
		double hmin=params.fget("hmin");
		double hmax=params.fget("hmax");
		int volume=params.iget("volume");
		double d=params.fget("dimension");
		
		for(double h=hmin;h<hmax;h+=0.0001){
			double ff=0;
			double lc=criticalLength(h);
			if(h>0 || lc>volume)lc=volume;
			clength.accum(h, lc);
			
			for(double l=1;l<lc;l+=0.001){
				ff+=Math.exp((-1/T)*(h*l+Math.pow(l,1-1/d)));
			}
			
			f.accum(h, ff);
		}
		
		double fe[]=f.copyData();
		double feprime[]=new double[fe.length];
		for(int i=0;i<fe.length/2-1;i++){
			feprime[2*i]=fe[2*i];
			feprime[2*i+1]=fe[2*i+3]-fe[2*i+1];
			fprime.accum(feprime[2*i],-feprime[2*i+1]);
			g.accum(-T*feprime[2*i+1], fe[2*i+1]);
		}

		Job.animate();
	}
	
	double criticalLength(double h){
		double d=2.0;
		return Math.pow((d-1)/(2*d*Math.abs(h)),d);		
	}

	public static void main(String[] args) {
		new Control(new cdmFreeEnergy(), "CDM Free Energy");
	}

}

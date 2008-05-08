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
	
	public void load(Control c){
		params.add("T",1.0);
		params.add("hmin", -0.1);
		params.add("hmax", 0.1);
		params.add("lmax", 1000);
		c.frame(F);
	}
	
	@Override
	public void animate() {
		F.registerPoints("Free Energy", f, Color.RED);
	}

	@Override
	public void clear() {
		f.clear();
		F.clear();
	}

	@Override
	public void run() {
		double T=params.fget("T");
		double hmin=params.fget("hmin");
		double hmax=params.fget("hmax");
		int LMAX=params.iget("lmax");
		for(double h=hmin;h<hmax;h+=0.0001){
			double ff=0;
			int lmax=0;
			if(h<0)lmax=(int)criticalLength(h);
			else lmax=LMAX;
			if(lmax>LMAX) lmax=LMAX;
			for(int l=1;l<lmax;l++){
				ff+=Math.exp((-1/T)*(h*l+Math.pow(l,0.5)));
			}
			f.accum(h, ff);
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

package rachele.OFC1d;

import java.awt.Color;

import kip.util.Random;
import scikit.dataset.Accumulator;
import scikit.graphics.dim2.Geom2D;
import scikit.graphics.dim2.Plot;
import scikit.jobs.Control;
import scikit.jobs.Job;
import scikit.jobs.Simulation;

public class OFC1dApp extends Simulation{


	double [] stress;
	boolean [] failedSites;
	boolean eqActive;
	int L;
	int failCount;
	int maxStressSite;
	double time;
	double alpha;
	
	Plot stressPlot = new Plot("Stress");
	
	protected Random random = new Random();
	
	
	public static void main(String[] args) {
		new Control(new OFC1dApp(), "OFC Damage Model Scaling Only");
	}
	
	public void animate() {
		params.set("Time",time);
		params.set("Fails", failCount);
		Accumulator acc = new Accumulator();
		Accumulator failAcc = new Accumulator();		
		for (int i=0; i<L; i++) acc.accum((double)(i), stress[i]);
		
		stressPlot.clearDrawables();

		stressPlot.addDrawable(
				Geom2D.line(0, 0.5, L, 0.5, Color.BLUE));
		stressPlot.addDrawable(
				Geom2D.line(0, 1.0, L, 1.0, Color.BLUE));
		for(int i = 0;i<L; i++){
			if(failedSites[i]){
				for(int j = 0; j < 20; j++){
					failAcc.accum((double)i+((double)j-10.0)/20.0, stress[i]);
				}

			}
		}
		if(eqActive){
			stressPlot.registerBars("Failed", failAcc, Color.red);
		}else{
			stressPlot.registerBars("Failed", failAcc, Color.green);		
		}
		stressPlot.registerBars("Stress", acc, Color.red);

	}

	public void clear() {
		
	}

	public void load(Control c) {
		c.frame(stressPlot);
		params.add("L", 16);
		params.add("alpha", 0.2);
		params.add("seed", 1);
		params.add("Fails");
		params.add("Time");
	}


	public void run() {
		L=params.iget("L");
		stress = new double [L];
		failedSites = new boolean [L];
		eqActive = false;
		alpha = params.fget("alpha");
		random.setSeed(params.iget("seed"));
		for (int i=0; i<L; i++)	stress[i] = random.nextDouble()*0.5+0.5;
		findMaxStressSite();

		time = 0;
		double stressStep = 0.01;
		
		
		while(true){
			if(eqActive){
				failCount +=1;
				failedSites[maxStressSite]=true;
				failSecondarySite();
				findMaxStressSite();
				if(stress[maxStressSite]>=1.0){
					eqActive = true;
				}else{
					eqActive = false;
				}
			}else{
				failCount = 0;
				double stressGap = 1.0-stress[maxStressSite];
				if(stressGap <= stressStep){
					for(int i = 0;i<L; i++){
						failedSites[i]=false;
					}
					failedSites[maxStressSite]=true;
					failCount +=1;
					for(int i = 0;i<L; i++){
						stress[i] += stressGap;
					}
					time += stressGap;
					stress[maxStressSite]=0.5;
					stress[(maxStressSite+1)%L]+=.25*(1-alpha);
					stress[(maxStressSite-1+L)%L]+=.25*(1-alpha);
					findMaxStressSite();
					if(stress[maxStressSite]>=1.0){
						eqActive = true;
					}else{
						eqActive = false;
					}
				}else{
					for(int i = 0;i<L; i++){
						stress[i] += stressStep;
					}
					time+=1.0;
				}
			}
		Job.animate();
		}

		
	}
	
	void findMaxStressSite(){
		maxStressSite = 0;
		for (int i=1; i<L; i++){
			if (stress[i]>stress[maxStressSite]) maxStressSite = i;
		}
	}
	void failSecondarySite(){
		double stressDrop = stress[maxStressSite]-0.5;
		stress[maxStressSite] -= stressDrop;
		stress[(maxStressSite+1)%L]+=(stressDrop/2.0)*(1-alpha);
		stress[(maxStressSite-1+L)%L]+=(stressDrop/2.0)*(1-alpha);
	}

}

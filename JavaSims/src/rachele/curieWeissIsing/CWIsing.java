package rachele.curieWeissIsing;

import kip.util.Random;
import scikit.jobs.params.Parameters;


/**
 * Curise-Weiss Ising model with Glauber dynamics.
 */
public class CWIsing {
	Random random = new Random();
	public int [] spin;
	public int L;
	int sumSpins;
	public double T, h, J;
	double dt, time;
	
	public CWIsing(Parameters params) {
		L = Integer.highestOneBit(params.iget("L"));
		params.set("L", L);
		spin = new int [L*L];
		random.setSeed(params.iget("Random seed"));
		setParameters(params);
	}
	
	public void setParameters(Parameters params) {
		
		dt = params.fget("dt");
		T  = params.fget("T");
		J  = params.fget("J", 1);
		h  = params.fget("h", 0.0);

	}
	
	public void randomizeField(double m) {
		if (m == 1 || m == -1) {
			for (int i = 0; i < L*L; i++)
				spin[i]=(int)m;
			sumSpins = (int)m*L*L;
		}
		else {
			int sum = 0;
			for (int i = 0; i < L*L; i++) {
				// p(s = +-1) = (1 +- m) / 2
				int s = (random.nextDouble() < (1+m)/2) ? 1 : -1;
				spin[i]=s;
				sum += s;
			}
			sumSpins = sum;
		}
	}
	
	public void step() {

		for (int cnt = 0; cnt < L*L*dt; cnt++) {
			//choose random spin
			int randomSpin =  (int)(random.nextDouble()*(double)(L*L));
			//calculate dE
			double dE = 2*spin[randomSpin]*(h + J*magnetization());
			//calculate prob of flipping for Glauber
			double prob = 1.0/(1.0+Math.exp(dE/T));
			//throw random number and choose flip
			if(random.nextDouble() < prob){
				spin[randomSpin] *= -1;
				sumSpins += 2*spin[randomSpin];
			}
		}
		time += dt;
	}
		
	
	public double magnetization() {
		return (double)sumSpins / (L*L);
	}
	
	public double time(){
		return time;
	}
	
	public double dTime(){
		return dt;
	}
	
	public int [] getField(int dx){
		return spin;
	}
	
	public void restartClock(){
		time = 0.0;
	}
	
}

package rachele.ising.dim2MC;

import scikit.jobs.params.Parameters;

/**
* Calculates energy metric for Ising model MC simulations.  The energy metric is defined as
* 
*  \Omega = (1/N)*\sum_{i=1}^N (ebar_i(t)-eave(t))
* 
*  where N = L*L, the size of the lattice, ebar is:
* 
*  ebar_i(t) = (1/t) \int_0^t dt' e_i(t')
* 
*  is the time averaged energy and 
*  
*  eave(t) = (1/N)*\sum_{i=1}^N(ebar_i(t)) is the averaged energy up to time t.
*  
* We allocate an ebar for each lattice point and calculate it at each
* time step.  Then we calc the averaged energy and metric at each time
* step. 
* 
*/
public class EnergyMetric {
	
	IsingLR sim;
	double dt, R;
	int L, N;
	double eAve; //averaged energy
	public double eMetric; 
	double [] eBar;  //time average of energy at one lattice point
	
	public EnergyMetric(IsingLR sim, Parameters params){
		this.sim = sim;
		this.L = sim.L;
		this.dt = sim.get_dt();
		this.R = sim.R;
		N = L*L;
		eBar = new double [N];
	}

	public void calculateMetric(){
		double ebarCalc;
		double pTime = sim.time() - dt;
		
		for (int i = 0; i < N; i++){
			ebarCalc = pTime*eBar[i];
			ebarCalc += sim.siteEnergy[i]*dt;
			eBar[i] = ebarCalc/sim.time();
		}
		
		double ebarSum = 0.0;
		for(int i = 0; i < N; i++) 
			ebarSum += eBar[i];
		eAve = ebarSum/N;
		
		double diffSum = 0;
		for(int i = 0; i < N; i++){
			diffSum += Math.pow(eBar[i] - eAve, 2);
		}
		eMetric = diffSum/N;
		
	}

}

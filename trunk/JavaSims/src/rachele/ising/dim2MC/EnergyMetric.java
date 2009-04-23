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
	double [] eBarP;
	int noParticles=0;
	double [] siteEnergy;
	double z;
	int rangeOffset;
	
	public EnergyMetric(IsingLR sim, Parameters params){
		this.sim = sim;
		this.L = sim.L;
		this.dt = sim.get_dt();
		this.R = sim.R;
		N = L*L;
		eBar = new double [N];
		siteEnergy = new double [L*L];
		z = (2*R+1)*(2*(R+rangeOffset)+1)-1;
		rangeOffset = sim.rangeOffset;
	}

	public void clearSums(){
		for (int i = 0; i < N; i++)
			eBar[i] = 0.0;
	}
	
	public void calculateMetric(){
		double ebarCalc;
		double pTime = sim.time() - dt;
		
		for (int i = 0; i < N; i++){
			ebarCalc = pTime*eBar[i];
			ebarCalc += siteEnergy[i]*dt;
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
	
	public void calculateParticleMetric(){
		double ebarCalc;
		double pTime = sim.time() - dt;
		
		for (int i = 0; i < noParticles; i++){
			ebarCalc = pTime*eBarP[i];
			ebarCalc += siteEnergy[sim.particleLoc[i]]*dt;
			eBarP[i] = ebarCalc/sim.time();
		}
		
		double ebarSum = 0.0;
		for(int i = 0; i < noParticles; i++)
			ebarSum += eBarP[i];
		eAve = ebarSum/noParticles;
		
		double diffSum = 0;
		for(int i = 0; i < noParticles; i++){
			diffSum += Math.pow(eBarP[i] - eAve, 2);
		}
		eMetric = diffSum/noParticles;		
	}
	
	public void calculateSiteEnergies(){
		for (int i = 0; i < N; i++){
			int x = i%L;
			int y = i/L;
			int s = sim.spins.get(x, y);
			siteEnergy[i] = -s*(sim.h + sim.J*(sim.sumWithOffset(x,y)-s)/z);
		}
	}
	
	public void calculateSiteEnergiesParticles(){
		for (int i = 0; i < noParticles; i ++){
			int x = i%L;
			int y = i/L;
			int s = sim.spins.get(x, y);
			siteEnergy[sim.particleLoc[i]] = -s*(sim.h + sim.J*(sim.sumWithOffset(x,y)-s)/z);
		}
	}
	
	public double calculateEnergyParticle(){
		double sum=0;
		for(int i = 0; i < noParticles; i++){
			sum += siteEnergy[sim.particleLoc[i]];
		}
		return sum;
	}
	
	public double calculateEnergy(){
		double sum=0;
		for(int i = 0; i < N; i++){
			sum += siteEnergy[i];
		}
		return sum;
	}
	
	/**
	 * Finds \Omega(t=0).  The time t=0 will not correspond to t=0 plotted
	 * above, but just needs to be taken at a time when the system is
	 * already in equilibrium (throw out the transient initial time.)
	 */
	public double findMetric0(){
		double eMetric0;
		double [] eBar0 = new double [N];
		
		for (int i = 0; i < N; i++)
			eBar0[i] = siteEnergy[i]*dt;
		
		double ebarSum = 0.0;
		for(int i = 0; i < N; i++) 
			ebarSum += eBar0[i];
		eAve = ebarSum/N;
		
		double diffSum = 0;
		for(int i = 0; i < N; i++){
			diffSum += Math.pow(eBar0[i] - eAve, 2);
		}
		eMetric0 = diffSum/N;
		return eMetric0;
//		System.out.println("Metric(t=0) = " + eMetric0);
		
	}

	public double findParticleMetric0(){
		double eMetric0;
		double [] eBar0 = new double [N];
		
		for (int i = 0; i < noParticles; i++)
			eBar0[i] = siteEnergy[sim.particleLoc[i]]*dt;
		
		double ebarSum = 0.0;
		for(int i = 0; i < noParticles; i++) 
			ebarSum += eBar0[i];
		eAve = ebarSum/noParticles;
		
		double diffSum = 0;
		for(int i = 0; i < noParticles; i++){
			diffSum += Math.pow(eBar0[i] - eAve, 2);
		}
		eMetric0 = diffSum/noParticles;
		return eMetric0;
//		System.out.println("Metric(t=0) = " + eMetric0);
		
	}
	
	public void initTrackEnergyMetric(int noP){
		noParticles = noP;
		eBarP = new double [noParticles];
	}

}

package chris.MD.TwoD;

import scikit.jobs.params.Parameters;
import chris.util.vector2d;

public class LennardJones extends InteractingSystem{
	
	static final double epsilon = 4; // this is really 4 epsilon, but we'll save on multiplications now
	static final double Vco = -epsilon*995904/244140625; // assumes a cutoff of 5/2 sigma
//	static final double Vco = 0;

	private double sigma, rco, rsm, rsd;

	public LennardJones(Parameters params){
		
		super(params);
		LJ_contructor(params);
	}
	
	public void LJ_contructor(Parameters params){
		
		sigma = params.fget("\u03C3");
		rco   = 2.5*sigma;
		rsm   = 2*sigma; 
		rsd   = (rco-rsm)*(rco-rsm)*(rco-rsm);
	}
	
	/*
	 * FORCES AND POTENTIALS DEPEND ON BC
	 * 
	 * 
	 */
	
	// ri are the "absolute" displacements of the two particles
	public vector2d force(vector2d r1, vector2d r2){
		
		return force(r1.x-r2.x, r1.y-r2.y);
	}
	
	// r is a relative displacement between two particles
	public vector2d force(vector2d r){
		
		return force(r.x, r.y);
	}
	
	// x and y are relative coordinate displacements between the two particles
	public vector2d force(double x, double y) {

		double r2 = x*x + y*y;
		if(r2 > rco*rco)
			return new vector2d();
				
		double sr2    = sigma*sigma/r2;
		double sr6    = sr2*sr2*sr2;
		double sr12   = sr6*sr6;
		double fmagXr = 6*epsilon*(2*sr12-sr6); // |r| \times |F_{LJ}| 

		if(r2 > rsm*rsm){
			double r = Math.sqrt(r2);
			fmagXr   = (fmagXr/r)*(1. - (r-rsm)*(r-rsm)*(3*rco - rsm - 2*r)/rsd) 
					   - (6*(r-rco)*(r-rsm)/rsd)*(epsilon*(sr12-sr6) - Vco);
			return new vector2d(x*fmagXr/r,y*fmagXr/r);  // |F_{LJ}| < x/r , y/r >
		}
		
		return new vector2d(x*fmagXr/r2,y*fmagXr/r2); // |F_{LJ}| < x/r , y/r >
	}

	public double potential(double r) {
		
		if(r > rco)
			return 0;
		
	    double sr   = sigma/r;
		double sr6  =  sr*sr*sr*sr*sr*sr;
		double sr12 =  sr6*sr6;
		
		if(r > rsm)
			return (1. - (r-rsm)*(r-rsm)*(3*rco - rsm - 2*r)/rsd)*(epsilon*(sr12-sr6) - Vco);
		
		return epsilon*(sr12-sr6) - Vco;
	}
	
	public double potential2(double r2) {
		
		if(r2 > rco*rco)
			return 0;
		
	    double sr2  = sigma*sigma/r2;
		double sr6  =  sr2*sr2*sr2;
		double sr12 =  sr6*sr6;
		
		if(r2 > rsm*rsm){
			double r = Math.sqrt(r2);
			return (1. - (r-rsm)*(r-rsm)*(3*rco - rsm - 2*r)/rsd)*(epsilon*(sr12-sr6) - Vco);

		}
		
		return epsilon*(sr12-sr6) - Vco;
	}

	

	
	
	
}

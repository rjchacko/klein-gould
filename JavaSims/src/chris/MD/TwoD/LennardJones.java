package chris.MD.TwoD;

import scikit.jobs.params.Parameters;
import chris.util.vector2d;

public class LennardJones extends InteractingSystem{
	

	static final double sigma   = 1.;
	static final double rco     = 2.5*sigma;
	static final double rsm     = 2*sigma; 
	static final double rsd     = (rco-rsm)*(rco-rsm)*(rco-rsm);
	static final double epsilon = 4; // this is really 4 epsilon
	static final double Vco     = epsilon*(Math.pow(sigma/rco,12)-Math.pow(sigma/rco,6));
	public LennardJones(Parameters params){
		
		super(params);
		LJ_contructor(params);
	}
	
	public void LJ_contructor(Parameters params){

		return;
	}
	
	
	// ri are the "absolute" displacements of the two particles
	public vector2d force(vector2d r1, vector2d r2){

		if(boundCond == BC.CLOSED)
			return force(r1.x-r2.x, r1.y-r2.y);

		double dx = r1.x-r2.x;
		double dy = r1.y-r2.y;
		if(dx*dx > 0.25*Lx*Lx)
			dx = -Math.signum(dx)*(Lx-Math.abs(dx));
		if(dy*dy > 0.25*Ly*Ly)
			dy = -Math.signum(dy)*(Ly-Math.abs(dy));
		return force(dx,dy);
	}
	
	// ri are the "absolute" displacements of the two particles
	public double potential(vector2d r1, vector2d r2){
		
		if(boundCond == BC.CLOSED)
			return potentialR2(vector2d.sub(r1, r2).length2());
		
		double dx = r1.x-r2.x;
		double dy = r1.y-r2.y;
		return potentialR2(Math.min(dx*dx,(Lx-dx)*(Lx-dx))+Math.min(dy*dy,(Ly-dy)*(Ly-dy)));
	}
	
	/*
	 * if this is the force of 1 ON 2 then m1 should be included as the coupling
	 * if it is the force of 2 ON 1 then m2 should be included as the coupling 
	 * maybe not . . . . ??
	 */
	// x and y are relative coordinate displacements between the two particles
	public vector2d force(double x, double y) {

//		double r2 = x*x + y*y;
//		if(r2 > rco*rco)
//			return new vector2d();
//				
//		double sr2    = sigma*sigma/r2;
//		double sr6    = sr2*sr2*sr2;
//		double sr12   = sr6*sr6;
//		double fmagXr = 6*epsilon*(2*sr12-sr6); // |r| \times |F_{LJ}| 
//
//		if(r2 > rsm*rsm){
//			double r = Math.sqrt(r2);
//			fmagXr   = (fmagXr/r)*(1. - (r-rsm)*(r-rsm)*(3*rco - rsm - 2*r)/rsd) 
//					   - (6*(r-rco)*(r-rsm)/rsd)*(epsilon*(sr12-sr6) - Vco);
//			return new vector2d(x*fmagXr/r,y*fmagXr/r);  // |F_{LJ}| < x/r , y/r >
//		}
//		
//		return new vector2d(x*fmagXr/r2,y*fmagXr/r2); // |F_{LJ}| < x/r , y/r >
		
		
		
		// no smoothing for debugging!
		
		double r2  = x*x + y*y;	
		if(r2 > rco*rco)
			return new vector2d();
		
		double sr2  = sigma*sigma/r2;
		double sr6  =  sr2*sr2*sr2;
		double sr12 = sr6*sr6;
		double fm = 6*epsilon*(2*sr12-sr6);
		return new vector2d(fm*x/r2,fm*y/r2);
	}

	// r2 is a relative displacement *squared* between two particles
	public double potentialR2(double r2) {
		
//		if(r2 > rco*rco)
//			return 0;
	
	    double sr2  = sigma*sigma/r2;
		double sr6  =  sr2*sr2*sr2;
		double sr12 =  sr6*sr6;
//		
//		if(r2 > rsm*rsm){
//			double r = Math.sqrt(r2);
//			return (1. - (r-rsm)*(r-rsm)*(3*rco - rsm - 2*r)/rsd)*(epsilon*(sr12-sr6) - Vco);
//
//		}
		// no smoothing  for debugging!
		return epsilon*(sr12-sr6) - Vco;
	}

	

	
	
	
}

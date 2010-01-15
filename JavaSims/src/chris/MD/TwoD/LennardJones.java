package chris.MD.TwoD;

import scikit.jobs.params.Parameters;
import chris.util.vector2d;

public class LennardJones extends InteractingSystem{
	
	static final double epsilon = 4; // this is really 4 epsilon, but we'll save on multiplications now
	static final double Vco = -epsilon*995904/244140625; // assumes a cutoff of 5/2 sigma
	
	private double sigma, rco;

	public LennardJones(Parameters params){
		
		super(params);
		LJ_contructor(params);
	}
	
	public void LJ_contructor(Parameters params){
		
		sigma = params.fget("\u03C3");
		rco   = 2.5*sigma;
		
	
		
	}
	
	public vector2d force(double x, double y) {

		double r2 = x*x + y*y;
		if(r2 > rco*rco)
			return new vector2d();
		
		double sr2    = sigma*sigma/r2;
		double sr6    = sr2*sr2*sr2;
		double sr12   = sr6*sr6;
		double fmagXr = 6*epsilon*(2*sr12-sr6); // |r| \times |F_{LJ}| 

		return new vector2d(x*fmagXr/r2,y*fmagXr/r2); // |F_{LJ}| < x/r , y/r >
	}

	public double potential(double r) {
		
		if(r > rco)
			return 0;
		
	    double sr   = sigma/r;
		double sr6  =  sr*sr*sr*sr*sr*sr;
		double sr12 =  sr6*sr6;
	
		return epsilon*(sr6-sr12) - Vco;
	}

	

	
	
	
}

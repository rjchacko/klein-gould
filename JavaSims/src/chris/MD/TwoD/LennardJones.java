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

	public vector2d force(vector2d deltaR) {
		
		double r2 = deltaR.length2();
		if(r2 > rco*rco)
			return new vector2d();
		
		double sr2  = sigma*sigma/r2;
		double sr6  =  sr2*sr2*sr2;
		double sr12 = sr6*sr6;
		double fm = 6*epsilon*(2*sr12-sr6);
		return new vector2d(fm*deltaR.x/r2,fm*deltaR.y/r2);
	}

	public double potential(vector2d deltaR) {
		
		double r2 = deltaR.length2();
		if(r2 > rco*rco)
			return 0;
		
	    double sr2  = sigma*sigma/r2;
		double sr6  =  sr2*sr2*sr2;
		double sr12 =  sr6*sr6;
		return epsilon*(sr12-sr6) - Vco;
	}

}

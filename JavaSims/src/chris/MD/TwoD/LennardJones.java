package chris.MD.TwoD;

import scikit.jobs.params.Parameters;
import chris.util.vector2d;

public class LennardJones extends InteractingSystem{
	
	static final double epsilon=1;

	public LennardJones(Parameters params){
		
		LJ_contructor(params);
	}
	
	public void LJ_contructor(Parameters params){
		
		
		
	}
	// this should be a system
	// it has particles, sigma, epsilon, temperature, etc
	
	
	
	public vector2d force(double r) {

		return null;
	}

	public double potential(double r) {

		return 0;
	}

	

	
	
	
}

package chris.MD.TwoD;

import scikit.jobs.params.Parameters;
import chris.util.vector2d;

public class IdealGas extends InteractingSystem{

	public IdealGas(Parameters params) {

		super(params);
	}

	public vector2d force(vector2d r1, vector2d r2) {
		
		return new vector2d();
	}

	public double potential(double r) {

		return 0;
	}

}

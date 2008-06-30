package chris.ofc;

import scikit.jobs.params.Parameters;

public class QuenchedDamage extends Damage2D {
	
	// Quench params
	
	private double SrQ0, SrQ[];
	private boolean SrQB;

	public QuenchedDamage(Parameters params) {
		super(params);
	
		Qconstructor(params);
		
	}

	public void Qconstructor(Parameters params){
		
		// set up Q noise
		if(params.fget("\u03C3_r Q-width") == 0){
			SrQ0 = 0;
			SrQB = false;
		}
		else{
			SrQ0 = params.fget("\u03C3_r Q-width");
			SrQB = true;
		}
		
		return;
	}
	
	public void foo(){
		System.out.println(SrQ0);
		if(SrQB){
			SrQ = null;
		}
		else{
			SrQ = new double[] {0.0 , 1.0 , 2.0 , 3.0};
		}
		System.out.println(SrQ);

		return;
	}
	
	
}

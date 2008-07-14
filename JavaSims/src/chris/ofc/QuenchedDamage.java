package chris.ofc;

import scikit.jobs.params.Parameters;

public class QuenchedDamage extends Damage2D {
	
	// Quench params
	
	private double SrQW, SrQ[], tSr0, tSrW;
	private boolean SrQB;
	
	
	public QuenchedDamage(Parameters params) {
		
		super(params);
		Qconstructor(params);
		
	}

	public void Qconstructor(Parameters params){
		
		// set up boolean and read-in parameters		
		SrQB = !((SrQW = params.fget("\u03C3_r Q-width")) == 0);
		tSr0 = getSr0();
		tSrW = getSrW();
		
		return;
	}
	
	public void Initialize(){
		
		super.Initialize();
		
		// set up SrQ
		SrQ = new double[getN()];
		
		if(SrQB){	// if quenched noise = true, set it up, else, do nothing
			for (int jj = 0 ; jj < getN() ; jj++){
				SrQ[jj] = tSr0 + SrQW*(rand.nextDouble() - 0.5);
			}
		}
		else{
			for (int jj = 0 ; jj < getN() ; jj++){
				SrQ[jj] = tSr0;
			}
		}
		
		return;
	}
	
	protected double nextSr(int site){
		
		return (SrQ[site] + tSrW*(rand.nextDouble()-0.5));
	}
	
	
}

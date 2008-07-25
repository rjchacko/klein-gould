package chris.ofc;

import chris.util.CopyArray;
import scikit.jobs.params.Parameters;

public class ForceSf extends Damage2D{

	private double SfNOW, dSfNOW;
	
	public ForceSf(Parameters params) {
	
		super(params);
		fConstructor(params);
	
	}

	public void fConstructor(Parameters params){
		
		dSfNOW = params.fget("d\u03C3_f");
		return;
	}

	public void Initialize(){
		
		super.Initialize();	
		
	}
	
	protected void InitArrays(){
		
		super.InitArrays();
			
		SfNOW = getSf0();
		
		// redistribute Sf and stress
		for (int jj = 0 ; jj < getN() ; jj++){
			Sf[jj] = getSr0() + (SfNOW - getSr0())*rand.nextDouble();
			setStress(jj, (1+10E-2)*getSr0());	// ad hoc!!!!
		}
		
		for (int jj = 0 ; jj < getN() ; jj++){
			Sf[jj] = getSr0() + (SfNOW - getSr0())*rand.nextDouble();
			setStress(jj, (1+10E-2)*getSr0());	// ad hoc!!!!
		}
		
	}
	
	protected void ForceFailure(){
			
		int[] temp = new int[getN()];
		int tempC = 0;
		
		if(getTime(-1) < 0){	// deal with equilibration
			seeds = null;
			return;		// SfNOW -= dSfNOW;	
		}
		
		SfNOW -= dSfNOW;	
		
		for (int jj = 0 ; jj < getN() ; jj++){
			Sf[jj] = getSr0() + (Sf[jj] - getSr0())*rand.nextDouble();
			if (getStress(jj) > Sf[jj]){
				temp[tempC++] = jj;
			}
		}
		
		seeds = CopyArray.copyArray(temp,tempC);
		
		ManageLives(seeds);
		setT1T2();	// for purposes of measuring the metric
		
		return;
	}
	
	public double getSfofT(){
		
		return SfNOW; 
	}
	
	
	
}

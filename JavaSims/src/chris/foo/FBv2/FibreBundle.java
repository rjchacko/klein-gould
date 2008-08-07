package chris.foo.FBv2;

import chris.foo.ofc.Damage2D;
import scikit.jobs.params.Parameters;

public class FibreBundle extends Damage2D{

	public FibreBundle(Parameters params) {
		super(params);
	}

	
	protected void InitArrays(){
		
		super.InitArrays();
		
		//double a = getSr0() + (getSf0()-getSr0())*0.001;	//can go no lower!!!!
		double b = getSf0();								// can go no higher!!!
		double a = 1.0*b;
			
		for (int jj = 0 ; jj < getN() ; jj++){
			
			//setStress(jj, getSr0() + (getSf0()-getSr0())*0.001);	// AD HOC!!
			setStress(jj,getSr0() + (getSf0() - getSr0())*rand.nextDouble());
			// create random number between [ 0 , b - a ] and then add a
			// [ 0 , b - a ] + a --> [a , b ]
			setSf(jj, a + (b-a)*rand.nextDouble());
			
		
		}
		
		
	}
	
}

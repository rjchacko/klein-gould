package chris.RandomWalker;

import scikit.numerics.fn.Function1D;

public class Pdf2step implements Function1D {

	private int dim;
	
	public Pdf2step(){
		
		Pdf2stepConstructor();
		return;
	}
	
	public void Pdf2stepConstructor(){
		
		dim = 2;
		return;
	}
	
	
	
	public double eval(double x) {
		
		if(dim==3){
			return x*x <= 1 ? 1 : 0;
		}
		else if(dim==1){
			return x*x == 1 ? 1 : 0;
		}
		else{
			return x*x <= 1 ? Math.pow(1-x*x,(dim-3)/2) : 0;
		}
	}

	public void setDim(int d){
		
		dim = d;
	}
	
	public int getDim(){
		
		return dim;
	}
	
}

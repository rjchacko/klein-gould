package rachele.ising.dim2;

import static scikit.numerics.Math2.j0;
import static scikit.numerics.Math2.j1;
import static scikit.numerics.Math2.jn;
import kip.util.Random;
import scikit.jobs.params.Parameters;

abstract public class AbstractIsing2D {
	public double L, R, T, dx, H;
	public int Lp;
	Random random = new Random();
	
	public static final double DENSITY = 1;
	
	// value of kR which minimizes potential
	public static final double KRcircle = 5.13562230184068255630140;
	public static final double KRsquare = 4.4934092;
	public static final double KRsquare2 = 11.6198;
	
	
	public double findVkCircle(double kR) {
		return (kR == 0 ? 1 : 2*j1(kR)/kR);
	}
	
	public double findVkSquare(double kRx, double kRy){
		double Vkx = (kRx == 0) ? 1 : Math.sin(kRx)/kRx;
		double Vky = (kRy == 0) ? 1 : Math.sin(kRy)/kRy;
		double Vk = Vkx*Vky;
		return Vk;
	}
	
	public double dpotential_dkR(double kR) {
		double kR2 = kR*kR;
		return (kR == 0) ? 0 : j0(kR)/kR - 2*j1(kR)/kR2  - jn(2,kR)/kR;
	}
	
	public int findKRSquareInt(double kRvalue){
		int kRint;
		int kRint1 = (int)((kRvalue*L)/(2*Math.PI*R));
		int kRint2 = kRint1+1;
		double kRValue1 = 2*Math.PI*kRint1*R/L;
		double kRValue2 = 2*Math.PI*kRint2*R/L;
		if(Math.abs(kRvalue-kRValue1) < Math.abs(kRvalue-kRValue2)){
			kRint = kRint1;
			System.out.println("kR square int= " + kRint + ", kRValue = " + kRValue1 + ", target= " + kRvalue + ", next guess = " + kRValue2);
		}else{
			kRint = kRint2;
			System.out.println("kR square int= " + kRint + ", kRValue = " + kRValue2 + ", target= " + kRvalue + ", next guess = " + kRValue1);
		}
		return kRint;
	}
	
	abstract public void readParams(Parameters params);
	abstract public void simulate();
	abstract public double[] coarseGrained();
	abstract public double time();
}

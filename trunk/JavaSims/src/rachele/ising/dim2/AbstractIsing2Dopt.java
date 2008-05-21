package rachele.ising.dim2;

import static scikit.numerics.Math2.j0;
import static scikit.numerics.Math2.j1;
import static scikit.numerics.Math2.jn;
import kip.util.Random;
import scikit.jobs.params.Parameters;

abstract public class AbstractIsing2Dopt {
	public double L, Rx, Ry, T, dx, H, J;
	public int Lp;
	Random random = new Random();
	
	public static final double DENSITY = 1;
	
	// value of kR which minimizes potential
	public static final double KRcircle = 5.13562230184068255630140;
	public static final double KRsquare = 4.4934092;
	
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
	
	abstract public void readParams(Parameters params);
	abstract public void simulate();
	abstract public double[] coarseGrained();
	abstract public double time();
}

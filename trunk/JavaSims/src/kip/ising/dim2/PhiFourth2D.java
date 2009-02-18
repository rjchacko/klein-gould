package kip.ising.dim2;

import static java.lang.Math.cos;
import static java.lang.Math.floor;
import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.Math.rint;
import static java.lang.Math.sqrt;
import static scikit.numerics.Math2.hypot;
import static scikit.numerics.Math2.j0;
import static scikit.numerics.Math2.j1;
import static scikit.numerics.Math2.jn;
import static scikit.numerics.Math2.sqr;
import static scikit.util.DoubleArray.copy;
import static scikit.util.DoubleArray.dot;
import static scikit.util.DoubleArray.variance;
import kip.util.Random;
import scikit.dataset.Accumulator;
import scikit.jobs.params.Parameters;
import scikit.numerics.fft.FFT2D;
import scikit.numerics.fn.Function2D;
import scikit.util.DoubleArray;


class PhiFourth2DConsts {
	// value of kR which minimizes j1(kR)/kR
	public static final double KR_SP = 5.13562230184068255630140;
	// S(k) ~ 1 / (V(kR_sp)/T+1)
	// => T_SP = - V(kR_sp) = - 2 j1(kR_sp) / kR_sp 
	public static final double T_SP = 0.132279487396100031736846;	
	
	public double potential(double kR) {
		return (kR == 0 ? 1 : 2*j1(kR)/kR);
	}
	
	public double dpotential_dkR(double kR) {
		double kR2 = kR*kR;
		return (kR == 0) ? 0 : j0(kR)/kR - 2*j1(kR)/kR2  - jn(2,kR)/kR;
	}
	
}

public class PhiFourth2D extends PhiFourth2DConsts {

	public double L, R, T, h, dx;
	Random random = new Random();
	
	int Lp;
	double t;
	double[] phi, phi_bar, del_phi;
	double[] background;
	double[] reactionVector;
	FFT2D fft;
	boolean noiselessDynamics = false;
	
	public double dt;
	public double Rx, Ry;
	public boolean rescaleClipped = false; // indicates saddle point invalid
	public double rms_dF_dphi;
	public double freeEnergyDensity;
	
	
	public PhiFourth2D(Parameters params) {
		random.setSeed(params.iget("Random seed", 0));
		
		Rx = Ry = params.fget("R");
		L = params.fget("L");
		T = params.fget("T");
		h = params.fget("h");
		dx = params.fget("dx");
		dt = params.fget("dt");
		noiselessDynamics = params.sget("Noise").equals("No");
		Lp = Integer.highestOneBit((int)rint(L/dx));
		dx = L / Lp;
		params.set("dx", dx);
		allocate();
		
		t = 0;
		for (int i = 0; i < Lp*Lp; i++)
			phi[i] = 0;
	}
	
	public double[] phi() {
		return phi;
	}
	
	public void halveResolution() {
		int old_Lp = Lp;
		double[] old_phi = phi; 
		Lp /= 2;
		dx *= 2.0;
		allocate();
		for (int y = 0; y < Lp; y++) {
			for (int x = 0; x < Lp; x++) {
				phi[y*Lp+x] = old_phi[2*y*old_Lp + 2*x];
				throw new Error("implement for bg");
			}
		}
	}
	
	public void doubleResolution() {
		int old_Lp = Lp;
		double[] old_phi = phi; 
		Lp *= 2;
		dx /= 2.0;
		allocate();
		for (int y = 0; y < Lp; y++) {
			for (int x = 0; x < Lp; x++) {
				phi[y*Lp+x] = old_phi[(y/2)*old_Lp + (x/2)];
				throw new Error("implement for bg");
			}
		}
	}
	
	public void readParams(Parameters params) {
		noiselessDynamics = params.sget("Noise").equals("No");
		T = params.fget("T");
		h = params.fget("h");
		dt = params.fget("dt");
	}
	
	
	public void initializeFieldWithRandomSeed() {
		copy(background, phi);
		
		for (int i = 0; i < Lp*Lp; i++) {
			double R = Rx;
			double x = dx*(i%Lp - Lp/2);
			double y = dx*(i/Lp - Lp/2);
			double r = sqrt(x*x+y*y);
			double mag = 0.8 / (1+sqr(r/R));
			phi[i] += mag*random.nextGaussian()/5;
		}
	}
	
	public void initializeFieldWithHexSeed() {
		copy(background, phi);
		
 		for (int i = 0; i < Lp*Lp; i++) {
			double R = Rx;
			double x = dx*(i%Lp - Lp/2);
			double y = dx*(i/Lp - Lp/2);
			double field = 0;
			double k = KR_SP/R;
			field = 0;
			field += cos(k * (1*x + 0*y));
			field += cos(k * (0.5*x + 0.5*sqrt(3)*y));
			field += cos(k * (-0.5*x + 0.5*sqrt(3)*y));

			double r = sqrt(x*x+y*y);
			double mag = 0.5 / (1+sqr(r/R));
			phi[i] += mag*field;
		}
	}
	
	public void useNoiselessDynamics(boolean b) {
		noiselessDynamics = b;
	}	
	
	public double phiVariance() {
		return variance(phi);
	}
		
	public void simulate() {
		fft.convolve(phi, phi_bar, new Function2D() {
			public double eval(double k1, double k2) {
				return potential(hypot(k1*Rx,k2*Ry));
			}
		});
		
		for (int i = 0; i < Lp*Lp; i++) {
			del_phi[i] = - dt*(phi_bar[i]+(T+T_SP)*phi[i] + phi[i]*phi[i]*phi[i] - h);
			if (!noiselessDynamics)
				del_phi[i] += sqrt(dt/(dx*dx))*noise();
		}
		
		rms_dF_dphi = 0;
		freeEnergyDensity = 0;
		for (int i = 0; i < Lp*Lp; i++) {
			rms_dF_dphi += sqr(del_phi[i] / dt);
			freeEnergyDensity += phi[i]*phi_bar[i]/2 + (T+T_SP)*phi[i]*phi[i]/2 + phi[i]*phi[i]*phi[i]*phi[i]/4 - h*phi[i];
			phi[i] += del_phi[i];
		}
		rms_dF_dphi = sqrt(rms_dF_dphi/(Lp*Lp));
		freeEnergyDensity /= (Lp*Lp);
		t += dt;
	}
	
	public void relaxInteraction() {
		Rx -= 0.1*dt*sqr(Rx)*dFdensity_dRx();
		Ry -= 0.1*dt*sqr(Ry)*dFdensity_dRy();
	}
	
	public void saddleStep() {
		double[] v = reactionVector;
		DoubleArray.sub(phi, background, v);
		
		double v2 = dot(v, v);
		double p1 = dot(v, phi);
		simulate();
		double p2 = dot(v, phi);
		double alpha = 0.5;
		for (int i = 0; i < Lp*Lp; i++) {
			phi[i] += (v[i]/v2) * (p1-p2) * (1+alpha);
			phi[i] = max(phi[i], -10);
			phi[i] = min(phi[i], +10);
		}
		
		iterateGrowthMode();
	}
	
	public void iterateGrowthMode() {
//		double[] v = growthEigenmode;
//		DoubleArray.normalize(v); // optional rescaling
//		
//		double[] v_bar = phi_bar; // use phi_bar as temporary space
//		fft.convolve(v, v_bar, new Function2D() {
//			public double eval(double k1, double k2) {
//				return potential(hypot(k1*Rx,k2*Ry));
//			}
//		});
//		double[] A_growth = phi_bar; // reuse the same temporary space
//		for (int i = 0; i < Lp*Lp; i++) {
//			A_growth[i] = v_bar[i]+(T/phi[i])*v[i];
//			v[i] += - dt*A_growth[i];
//		}
//		
//		growthEigenvalue = dot(v, A_growth) / dot(v, v);		
//		double[] zero = A_growth;
//		for (int i = 0; i < Lp*Lp; i++)
//			zero[i] = A_growth[i] - growthEigenvalue*v[i];
//		double mag2 = sqr(growthEigenvalue)*dot(v, v);
//		rms_growthEigenmode = sqrt(dot(zero, zero) / mag2);
	}
	
	public double dFdensity_dRx() {
		double[] dphibar_dR = phi_bar;
		fft.convolve(phi, phi_bar, new Function2D() {
			public double eval(double k1, double k2) {
				double kR = hypot(k1*Rx, k2*Ry);
				double dkR_dRx = k1 == 0 ? 0 : (k1*k1*Rx / kR);
				return dpotential_dkR(kR)*dkR_dRx;
			}
		});
		return DoubleArray.dot(phi, dphibar_dR) / (2*Lp*Lp);
	}
	
	public double dFdensity_dRy() {
		double[] dphibar_dR = phi_bar;
		fft.convolve(phi, phi_bar, new Function2D() {
			public double eval(double k1, double k2) {
				double kR = hypot(k1*Rx, k2*Ry);
				double dkR_dRy = k2 == 0 ? 0 : (k2*k2*Ry / kR);
				return dpotential_dkR(kR)*dkR_dRy;
			}
		});
		return DoubleArray.dot(phi, dphibar_dR) / (2*Lp*Lp);
	}

	public Accumulator newStructureAccumulator(double binWidth) {
		// round binwidth down so that it divides KR_SP without remainder.
		binWidth = KR_SP / floor(KR_SP/binWidth);
		return new Accumulator(binWidth);
	}
	
	public void accumulateStructure(final Accumulator sf) {
		fft.transform(phi, new FFT2D.MapFn() {
			public void apply(double k1, double k2, double re, double im) {
				double kR = hypot(k1*Rx, k2*Ry);
				if (kR > 0 && kR <= 4*KR_SP)
					sf.accum(kR, (re*re+im*im)/(L*L));
			}
		});
	}	
	
	public void saveFieldAsBackground() {
		copy(phi, background);
	}
	
	public int numColumns() {
		return Lp;
	}
	
	public double time() {
		return t;
	}
	
	
	private void allocate() {
		phi = new double[Lp*Lp];
		phi_bar = new double[Lp*Lp];
		del_phi = new double[Lp*Lp];
		background = new double[Lp*Lp];
		reactionVector = new double[Lp*Lp];
		
		fft = new FFT2D(Lp, Lp);
		fft.setLengths(L, L);
	}
	
	private double noise() {
		return noiselessDynamics ? 0 : random.nextGaussian();
	}

}

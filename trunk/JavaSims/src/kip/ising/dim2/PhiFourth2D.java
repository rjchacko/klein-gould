package kip.ising.dim2;

import static java.lang.Math.floor;
import static java.lang.Math.rint;
import static java.lang.Math.sqrt;
import static scikit.numerics.Math2.hypot;
import static scikit.numerics.Math2.j1;
import static scikit.numerics.Math2.sqr;
import kip.util.Random;
import scikit.dataset.Accumulator;
import scikit.jobs.params.Parameters;
import scikit.numerics.fft.FFT2D;
import scikit.numerics.fn.Function2D;


class PhiFourth2DConsts {
	// value of kR which minimizes j1(kR)/kR
	public static final double KR_SP = 5.13562230184068255630140;
	// S(k) ~ 1 / (V(kR_sp)/T+1)
	// => T_SP = - V(kR_sp) = - 2 j1(kR_sp) / kR_sp 
	public static final double T_SP = 0.132279487396100031736846;	
	
	Function2D potential1 = new Function2D() {
		public double eval(double x, double y) {
			return - (x*x + y*y);
		}
	};

	Function2D potential2 = new Function2D() {
		public double eval(double kx, double ky) {
			double k = hypot(kx, ky);
			return (k == 0 ? 1 : 2*j1(k)/k);
		}
	};
}

public class PhiFourth2D extends PhiFourth2DConsts {

	public double L, R, T, h, dx;
	Random random = new Random();
	
	int Lp;
	double t;
	double[] phi, phi_bar, del_phi;
	FFT2D fft;
	boolean noiselessDynamics = false;
	
	public double dt;
	public double rms_dF_dphi;
	public double freeEnergyDensity;
	
	
	public PhiFourth2D(Parameters params) {
		random.setSeed(params.iget("Random seed", 0));
		
		R = params.fget("R");
		L = R*params.fget("L/R");
		T = params.fget("T");
		h = params.fget("h");
		dx = R*params.fget("dx/R");
		dt = params.fget("dt");
		noiselessDynamics = params.sget("Noise").equals("No");
		Lp = Integer.highestOneBit((int)rint(L/dx));
		dx = L / Lp;
		params.set("dx/R", dx/R);
		allocate();
		
		t = 0;
		for (int i = 0; i < Lp*Lp; i++)
			phi[i] = 0;
	}
	
	public double[] phi() {
		return phi;
	}
	
	public void readParams(Parameters params) {
		noiselessDynamics = params.sget("Noise").equals("No");
		T = params.fget("T");
		h = params.fget("h");
		dt = params.fget("dt");
	}
	
	
	public void randomize() {
		for (int i = 0; i < Lp*Lp; i++) {
			phi[i] += sqrt(1/(R*R*dx*dx))*random.nextGaussian();
		}
	}
	
	public void useNoiselessDynamics(boolean b) {
		noiselessDynamics = b;
	}	
	
	public void simulate() {
		fft.convolve(phi, phi_bar, potential2);
		
		for (int i = 0; i < Lp*Lp; i++) {
			del_phi[i] = - dt*(phi_bar[i]+T*phi[i] + phi[i]*phi[i]*phi[i] - h);
			if (!noiselessDynamics)
				del_phi[i] += sqrt(dt/(R*R*dx*dx))*noise();
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

	public Accumulator newStructureAccumulator(double binWidth) {
		// round binwidth down so that it divides KR_SP without remainder.
		binWidth = KR_SP / floor(KR_SP/binWidth);
		return new Accumulator(binWidth);
	}
	
	public void accumulateStructure(final Accumulator sf) {
		fft.transform(phi, new FFT2D.MapFn() {
			public void apply(double k1, double k2, double re, double im) {
				double kmag = hypot(k1, k2);
				if (kmag > 0 && kmag <= 4)
					sf.accum(kmag, (re*re+im*im)/(L*L));
			}
		});
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
		
		fft = new FFT2D(Lp, Lp);
		fft.setLengths(L, L);
	}
	
	private double noise() {
		return noiselessDynamics ? 0 : random.nextGaussian();
	}

}

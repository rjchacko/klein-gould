package kip.ising.dim2;

import static java.lang.Math.rint;
import static java.lang.Math.sqrt;
import static scikit.numerics.Math2.hypot;
import static scikit.numerics.Math2.sqr;
import kip.util.Random;
import scikit.dataset.Accumulator;
import scikit.jobs.params.Parameters;
import scikit.numerics.fft.FFT2D;
import scikit.numerics.fn.Function2D;


class PhiFourth2DConsts {
	Function2D potential1 = new Function2D() {
		public double eval(double kx, double ky) {
			return - (kx*kx + ky*ky);
		}
	};
}

public class PhiFourth2D extends PhiFourth2DConsts {

	public double L, R, T, h, dx;
	Random random = new Random();
	
	int Lp;
	double t;
	double[] phi, laplace_phi, del_phi;
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
			phi[i] = sqrt(1/(dx*dx))*random.nextGaussian();
			phi[i] += 1;
		}
		
		System.out.println(Lp);
		for (int x = 0; x < Lp; x++) {
			for (int y = 0; y < Lp; y++) {
				int i = y*Lp + x;
				phi[i] = Math.sin(0.5*(x+32)); // hypot(x*dx, y*dx); 
			}
		}
	}
	
	public void useNoiselessDynamics(boolean b) {
		noiselessDynamics = b;
	}	
	
	private void laplaceOperator(double[] phi, double[] phi_bar) {
		for (int x = 0; x < Lp; x++) {
			for (int y = 0; y < Lp; y++) {
				int xp = (x+1)%Lp;
				int xm = (x-1+Lp)%Lp;
				int yp = (y+1)%Lp;
				int ym = (y-1+Lp)%Lp;
				phi_bar[y*Lp+x] = (-4*phi[y*Lp+x] + phi[yp*Lp+x] + phi[ym*Lp+x] + 
									phi[y*Lp+xp] + phi[y*Lp+xm]) / (dx*dx); 
			}
		}
	}
	
	public void simulate() {
		fft.convolve(phi, laplace_phi, potential1);
//		laplaceOperator(phi, laplace_phi);
		
		for (int i = 0; i < Lp*Lp; i++) {
//			del_phi[i] = - dt*(phi_bar[i]+T*phi[i] + phi[i]*phi[i]*phi[i] - h);
			del_phi[i] = - dt*(-R*R*laplace_phi[i] + (T+sqr(phi[i]))*phi[i] - h);
			if (!noiselessDynamics)
				del_phi[i] += sqrt(dt/(dx*dx))*noise();
		}
		
//		rms_dF_dphi = 0;
//		freeEnergyDensity = 0;
		for (int i = 0; i < Lp*Lp; i++) {
//			rms_dF_dphi += sqr(del_phi[i] / dt);
			freeEnergyDensity += -phi[i]*(R*R*laplace_phi[i])/2 + sqr(T+sqr(phi[i]))/4 - h*phi[i];
			phi[i] += del_phi[i];
		}
//		rms_dF_dphi = sqrt(rms_dF_dphi/(Lp*Lp));
//		freeEnergyDensity /= (Lp*Lp);
		t += dt;
	}

	public Accumulator newStructureAccumulator(double binWidth) {
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
		laplace_phi = new double[Lp*Lp];
		del_phi = new double[Lp*Lp];
		
		fft = new FFT2D(Lp, Lp);
		fft.setLengths(L, L);
	}
	
	private double noise() {
		return noiselessDynamics ? 0 : random.nextGaussian();
	}

}

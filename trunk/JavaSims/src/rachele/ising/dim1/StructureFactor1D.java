package rachele.ising.dim1;

import static java.lang.Math.PI;
//import static java.lang.Math.sin;
//import static java.lang.Math.sqrt;
import scikit.dataset.Accumulator;
//import scikit.numerics.fft.ComplexDouble2DFFT;
import scikit.numerics.fft.managed.ComplexDoubleFFT;
import scikit.numerics.fft.managed.ComplexDoubleFFT_Mixed;

public class StructureFactor1D {
	//RealDoubleFFT_Radix2 fft;	// Object to perform transforms

	ComplexDoubleFFT fft;
	double[] fftData;       // Fourier transform data
	int Lp;                 // # elements per side
	double L;               // the actual system length, L = Lp*dx, where dx is lattice spacing
	double R;               // characteristic length.  x-axis is k*R.
	double kRmin, kRmax;
	Accumulator acc;

	public StructureFactor1D(int Lp, double L, double R, double kRbinWidth) {
		this.Lp = Lp;
		this.L = L;
		this.R = R;
		
		kRmin = (2*PI*2/L)*R; // explicitly exclude constant (k=0) mode
		kRmax = (2*PI*(Lp/2)/L)*R;
		acc = new Accumulator(kRbinWidth);
		//fft = new RealDoubleFFT_Radix2(Lp);
		fft = new ComplexDoubleFFT_Mixed(Lp);
		fftData = new double[2*Lp];
	}
	
	public Accumulator getAccumulator() {
		return acc;
	}
	
	public void accumulate(double[] xs) {
		for (int i = 0; i < Lp; i++) {
			fftData[2*i+0] = xs[i];
			fftData[2*i+1] = 0;
		}
		fft.transform(fftData);
		for (int x = -Lp/2; x < Lp/2; x++) {
			double kR = (2*PI*x/L) * R;
			if (kR >= kRmin && kR <= kRmax) {
				int i = (x + Lp) % Lp;
				double re = fftData[2*i];
				double im = fftData[2*i+1];
				acc.accum(kR, (re*re + im*im)/(L));
				}
		}
	}
	
}

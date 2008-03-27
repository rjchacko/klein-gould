package rachele.ising.dim1;

import kip.util.Random;
import scikit.dataset.Accumulator;
import scikit.jobs.params.Parameters;
import scikit.numerics.fft.managed.ComplexDoubleFFT;
import scikit.numerics.fft.managed.ComplexDoubleFFT_Mixed;
import static java.lang.Math.PI;
import static java.lang.Math.log;
import static java.lang.Math.rint;
import static java.lang.Math.sin;
import static java.lang.Math.sqrt;
import static scikit.numerics.Math2.atanh;
import static scikit.util.DoubleArray.*;

//import java.util.*;

public class FieldIsing1D{
	public int Lp;
	public double dt, t;
	public double[] phi, F;
	double DENSITY;
	double [] phi_bar, del_phi;
	boolean modelA, noise;
	
	public double L, R, T, J, dx, H;
	Random random = new Random();
	
	//public static final double DENSITY = -0;
	public static final double KR_SP = 4.4934102;
	
	ComplexDoubleFFT fft;
	private double[] fftScratch;
	public double freeEnergyDensity;
	
	Accumulator freeEngAcc;
	
	public FieldIsing1D(Parameters params) {
		random.setSeed(params.iget("Random seed", 0));
		
		R = params.fget("R");
		L = R*params.fget("L/R");
		T = params.fget("T");
		dx = R/params.fget("R/dx");
		dt = params.fget("dt");
		J = params.fget("J");
		DENSITY = params.fget("Density");
		H = params.fget("H");
		Lp = Integer.highestOneBit((int)rint((L/dx)));
		dx = L / Lp;
		double RoverDx = R/dx;
		params.set("R/dx", RoverDx);
		params.set("Lp", Lp);
		if (params.sget("Model").equals("A"))
			modelA = true;
		else
			modelA = false;
		t = 0;

		phi = new double[Lp];
		phi_bar = new double[Lp];
		del_phi = new double[Lp];
		F = new double [Lp];
		
		fftScratch = new double[2*Lp];
		fft = new ComplexDoubleFFT_Mixed(Lp);
		freeEngAcc = new Accumulator(dt);
		
		for (int i = 0; i < Lp; i++)
			phi[i] = DENSITY;
	}
	
	public void readParams(Parameters params) {
		T = params.fget("T");
		J = params.fget("J");
		dt = params.fget("dt");
		R = params.fget("R");
		H = params.fget("H");
		L = R*params.fget("L/R");
		dx = R/params.fget("R/dx");
		Lp = Integer.highestOneBit((int)rint((L/dx)));
		dx = L / Lp;
		params.set("R/dx", R/dx);
		if (params.sget("Model").equals("A"))
			modelA = true;
		else
			modelA = false;
		params.set("DENSITY", mean(phi));
		if (params.sget("Noise").equals("On"))
			noise = true;
		else
			noise = false;		
	}
	
	public double time() {
		return t;
	}
	
	public Accumulator getFreeEngAcc() {
		return freeEngAcc;
	}
	
	public void measureFreeEng(){
		convolveWithRange(phi, phi_bar, R);
		freeEnergyDensity = 0;
		for (int i = 0; i < Lp; i ++){
			double potential = (phi[i]*phi_bar[i])/2.0;
			double entropy = -((1.0 + phi[i])*log(1.0 + phi[i]) +(1.0 - phi[i])*log(1.0 - phi[i]))/2.0;
			F[i] = potential - H*phi[i] - T*entropy; 
			freeEnergyDensity += F[i];
		}
		freeEnergyDensity /= (double)Lp;
		freeEngAcc.accum(t, freeEnergyDensity);
	}
	
	void convolveWithRange(double[] src, double[] dest, double R) {
		// write real and imaginary components into scratch
		for (int i = 0; i < Lp; i++) {
			fftScratch[2*i+0] = src[i];
			fftScratch[2*i+1] = 0;
		}
		
		// multiply real and imaginary components by the fourier transform
		// of the potential, V(k), a real quantity.  this corresponds to
		// a convolution in "x" space.
		fft.transform(fftScratch);
		for (int x = -Lp/2; x < Lp/2; x++) {
			double kR = (2*PI*x/L) * R;
			int i = (x + Lp) % Lp;
			double V = (kR == 0 ? 1 : sin(kR)/kR);
			fftScratch[2*i+0] *= J*V;
			fftScratch[2*i+1] *= J*V;
		}
		fft.backtransform(fftScratch);
		
		// after reverse fourier transformation, return the real result.  the
		// imaginary component will be zero.
		for (int i = 0; i < Lp; i++) {
			dest[i] = fftScratch[2*i+0] / Lp;
		}
	}
	
	public double[] copyField() {
		double ret[] = new double[Lp];
		for (int i = 0; i < Lp; i++)
			ret[i] = phi[i];
		return ret;
	}
	
	public void simulate() {
		convolveWithRange(phi, phi_bar, R);

		if (modelA){
			for (int i = 0; i < Lp; i++) {
				//del_phi[i] = - dt*( phi_bar[i]-H-T*log(1.0-phi[i])/2.0+T*log(1.0+phi[i])/2.0) + sqrt(dt*2*T/dx)*random.nextGaussian();
				del_phi[i] = - dt*(phi_bar[i]-H + T*atanh(phi[i]));
				if(noise)
					del_phi[i] += sqrt(dt*2*T/dx)*random.nextGaussian();
			}
			//double mu = mean(del_phi)-(DENSITY-mean(phi));
			for (int i = 0; i < Lp; i++) {
				phi[i] += del_phi[i];	
			}		
		}else{
			for (int i = 0; i < Lp; i++) {
				del_phi[i] = - dt*( phi_bar[i]-H-T*log(1.0-phi[i])+T*log(1.0+phi[i]));
				if(noise)
					del_phi[i] += sqrt(dt*2*T/dx)*random.nextGaussian();
			}
			double mu = mean(del_phi)-(DENSITY-mean(phi));
			for (int i = 0; i < Lp; i++) {
				phi[i] += del_phi[i] - mu;	
			}			
		}
		measureFreeEng();
		t += dt;
	}
}

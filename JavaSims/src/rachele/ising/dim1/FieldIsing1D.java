package rachele.ising.dim1;

import kip.util.Random;
//import rachele.util.FourierTransformer;
import scikit.dataset.Accumulator;
import scikit.jobs.params.Parameters;
import scikit.numerics.fft.managed.ComplexDoubleFFT;
import scikit.numerics.fft.managed.ComplexDoubleFFT_Mixed;
import static java.lang.Math.PI;
import static java.lang.Math.log;
import static java.lang.Math.pow;
import static java.lang.Math.rint;
import static java.lang.Math.sin;
import static java.lang.Math.sqrt;
import static scikit.numerics.Math2.atanh;
import static scikit.numerics.Math2.j1;
import static scikit.numerics.Math2.sqr;
import static scikit.util.DoubleArray.*;

/**
* This is for only Model A right now 
*/
public class FieldIsing1D{
	public int Lp;
	public double dt = 1;
	public double t, noiseParam;
	public double[] phi, F;
	double DENSITY;
	double [] phi_bar, del_phi;
//	boolean noise;
	
	public double L, R, T, J, dx, H, ampFactor;
	Random random = new Random();
	
	public static final double KR_SP = 4.4934102;
	
	ComplexDoubleFFT fft;
	private double[] fftScratch;
	public double freeEnergyDensity;
	
	Accumulator freeEngAcc;
	
	public FieldIsing1D(Parameters params) {
		random.setSeed(params.iget("Random Seed"));
		
		R = params.fget("R");
		L = R*params.fget("L/R");
		T = params.fget("T");
		dx = R/params.fget("R/dx");
		J = params.fget("J");
		DENSITY = params.fget("Density");
		H = params.fget("H");
		Lp = Integer.highestOneBit((int)rint((L/dx)));
		dx = L / Lp;
		double RoverDx = R/dx;
		params.set("R/dx", RoverDx);
		params.set("Lp", Lp);
		t = 0;
		noiseParam = params.fget("Noise");
		
		phi = new double[Lp];
		phi_bar = new double[Lp];
		del_phi = new double[Lp];
		F = new double [Lp];
		
		fftScratch = new double[2*Lp];
		fft = new ComplexDoubleFFT_Mixed(Lp);
		freeEngAcc = new Accumulator(dt);
		
		for (int i = 0; i < Lp; i++)
			phi[i] = DENSITY+ noiseParam*random.nextGaussian()*sqrt((1-DENSITY*DENSITY)/(dx*dx));;
	}
	
	public void readParams(Parameters params) {
		T = params.fget("T");
		J = params.fget("J");
		R = params.fget("R");
		H = params.fget("H");
		L = R*params.fget("L/R");
		dx = R/params.fget("R/dx");
		Lp = Integer.highestOneBit((int)rint((L/dx)));
		dx = L / Lp;
		params.set("R/dx", R/dx);
		params.set("DENSITY", mean(phi));
		noiseParam = params.fget("Noise");
//		if (params.sget("Noise").equals("On"))
//			noise = true;
//		else
//			noise = false;		
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
		for (int i = 0; i < Lp; i++) {
			double Lambda = sqr(1 - phi[i]*phi[i]);	
			del_phi[i] = - dt*Lambda*(phi_bar[i]-H+ T*atanh(phi[i]));
			del_phi[i] += sqrt(Lambda*dt*2*T/dx)*random.nextGaussian()*noiseParam;
		}
		for (int i = 0; i < Lp; i++)
			phi[i] += del_phi[i];	
		measureFreeEng();
		t += dt;
	}

	public void simulateGlauber() {
		
		convolveWithRange(phi, phi_bar, R);	
		for (int i = 0; i < Lp; i++){
			double arg = (phi_bar[i]/T + H/T);
			double tanh = Math.tanh(arg);
			double driftTerm = dt*(tanh-phi[i]);
			double noisePre = sqrt(2-pow(Math.tanh(arg),2)-pow(phi[i],2));
			double noiseTerm = noiseParam*noisePre*random.nextGaussian()*sqrt(dt*2/(dx*dx));
			del_phi[i] = driftTerm + noiseTerm;
		}

		for (int i = 0; i < Lp*Lp; i++) 
			phi[i] += del_phi[i];			
		t += dt;	
	}
	
	/**
	 * This method exists to find a stripe
	 * solution for the 2D circle interaction,
	 * which is not the same as the 1D solution
	 * with a step function potential.  To account for the 
	 * the "lost square corners" we have an effective
	 * potential that decreases at the edges of the 
	 * interaction range:  the interaction strength goes
	 * like sqrt(R*R-r*r).  I can't find a FT for this
	 * (can't do the integral) so I have to do it 
	 * numerically.
	 */
	public void simulateCircle(double ampFactor){
		//convolveCircleInteraction(phi, circleIntFT, phi_bar, R);
		//FourierTransformer FT = new FourierTransformer(Lp);
		//phi_bar = FT.convolve1D(phi, circleIntFT);
		//System.out.println("sim circle");
		convolveCircleInteraction(phi, phi_bar, R);
		for (int i = 0; i < Lp; i++) {
			double Lambda = sqr(1 - phi[i]*phi[i]);	
			del_phi[i] = - dt*Lambda*(ampFactor*phi_bar[i]-H + T*atanh(phi[i]));
			del_phi[i] += noiseParam*sqrt(Lambda*dt*2*T/dx)*random.nextGaussian();
		}
		for (int i = 0; i < Lp; i++)
			phi[i] += del_phi[i];	
		measureFreeEng();
		t += dt;
	}
	
	void convolveCircleInteraction(double[] src, double[] dest, double R) {
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
			double V = (kR == 0 ? PI/2 : PI*j1(kR)/kR);
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
	
}

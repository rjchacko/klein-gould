package rachele.ising.dim1;

import kip.util.Random;
//import rachele.util.FourierTransformer;
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
	public double dt;
	public double t, noiseParam;
	public double[] phi, F, phi_k;
	double DENSITY;
	double [] phi_bar, del_phi;
	public double L, R, T, J, dx, H, ampFactor;
	Random random = new Random();
	
	public static final double KR_SP = 4.4934102;
	
	ComplexDoubleFFT fft;
	private double[] fftScratch;
	public double freeEnergyDensity;

	public FieldIsing1D(Parameters params) {
		random.setSeed(params.iget("Random Seed"));
		
		R = params.fget("R");
		L = R*params.fget("L/R");
		T = params.fget("T");
		dx = R/params.fget("R/dx");
		J = params.fget("J");
		DENSITY = params.fget("Density");
		H = params.fget("H");
		dt = params.fget("dt");
		Lp = Integer.highestOneBit((int)rint((L/dx)));
		dx = L / Lp;
		double RoverDx = R/dx;
		params.set("R/dx", RoverDx);
		params.set("Lp", Lp);
		t = 0;
		noiseParam = params.fget("Noise");
		
		phi = new double[Lp];
//		phi_k = new double [Lp*2];
		phi_bar = new double[Lp];
		del_phi = new double[Lp];
		F = new double [Lp];
		
		fftScratch = new double[2*Lp];
		fft = new ComplexDoubleFFT_Mixed(Lp);
		
		for (int i = 0; i < Lp; i++)
			phi[i] = DENSITY+ noiseParam*random.nextGaussian()*sqrt((1-DENSITY*DENSITY)/(dx*dx));
		phi_k = transform(phi);
	}
	
	public void readParams(Parameters params) {
		T = params.fget("T");
		J = params.fget("J");
		R = params.fget("R");
		H = params.fget("H");
		L = R*params.fget("L/R");
		dx = R/params.fget("R/dx");
		dt = params.fget("dt");
		Lp = Integer.highestOneBit((int)rint((L/dx)));
		dx = L / Lp;
		params.set("R/dx", R/dx);
		params.set("DENSITY", mean(phi));
		noiseParam = params.fget("Noise");
		
	}
	
	public double time() {
		return t;
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
	
	public void simulateConserved(){
		convolveWithRange(phi, phi_bar, R);
		double drift [] = new double [Lp];
		for (int i = 0; i < Lp; i++) {
			double dF_dPhi = 0;
			dF_dPhi = (phi_bar[i] + T*scikit.numerics.Math2.atanh(phi[i]));
			drift[i] = - dt*dF_dPhi;
		}

		double [] k2ft_drift = transform(drift);
		for (int i = 0; i < Lp; i++) {
			double kValue = (2*PI*i*R/L);
			k2ft_drift[2*i] *= (kValue*kValue);
			k2ft_drift[2*i+1] *= (kValue*kValue);
		}
		del_phi = backtransform(k2ft_drift);
		
		for (int i = 0; i < Lp; i++) {
			phi[i] += del_phi[i] + sqrt(dt*2*T/dx)*random.nextGaussian()*noiseParam;
		}	
		
		t += dt;

	}
	
	public void simulateConservedWithMobility(){
		// dphi/dt = grad dot mobility * grad(del F/ del phi)
		// impliment as IFT(k*FT(mobility *IFT(k(FT(del F / del phi)))))
		// or
		// k*f_(mob*{if_[k(f_dfdp)]})
		convolveWithRange(phi, phi_bar, R);
		double drift [] = new double [Lp];
		for (int i = 0; i < Lp; i++) {
			double dF_dPhi = 0;
			dF_dPhi = (phi_bar[i] + T*scikit.numerics.Math2.atanh(phi[i]));
			drift[i] = - dt*dF_dPhi;
		}
		double [] kft_drift = transform(drift);
		for (int i = 0; i < Lp; i++) {
			double kValue = (2*PI*i*R/L);
			kft_drift[2*i] *= (kValue);
			kft_drift[2*i+1] *= (kValue);
		}
		double [] m_kdrift = backtransform(kft_drift);
		for (int i = 0; i < Lp; i++) 
			m_kdrift[i] *= Math.pow(1-phi[i]*phi[i],2);
		double [] k_mkdrift = transform(m_kdrift);
		for (int i = 0; i < Lp; i++) {
			double kValue = (2*PI*i*R/L);
			k_mkdrift[2*i] *= (kValue);
			k_mkdrift[2*i+1] *= (kValue);
		}		
		del_phi = backtransform(k_mkdrift);
		for (int i = 0; i < Lp; i++) {
			phi[i] += del_phi[i];// + sqrt(dt*2*T/dx)*random.nextGaussian()*noiseParam;
		}	
		t += dt;
	}
	
	public void simulateConseveredSemiImp(){
//		double [] atanh = new double [Lp];
//		for (int i = 0; i < Lp; i++){
//			atanh[i] = T*scikit.numerics.Math2.atanh(phi[i]);
//		}
//		double [] atanh_k = transform(atanh);
//		for (int x = -Lp/2; x < Lp/2; x++) {
//			double kR = (2*PI*x*R/L);
//			int i = (x + Lp) % Lp;
//			double V = (kR == 0 ? 1 : sin(kR)/(kR));
//			double k = (2*PI*x/L);
//			phi_k[2*i] = (phi_k[2*i] + dt*k*k*R*R*T*atanh_k[2*i])/(1.0-k*k*R*R*V*dt); 
//			phi_k[2*i+1] = (phi_k[2*i+1] + dt*k*k*R*R*T*atanh_k[2*i+1])/(1.0-k*k*R*R*V*dt); 
//		}
//		phi = backtransform(phi_k);
		
		convolveWithRange(phi, phi_bar, R);
		double drift [] = new double [Lp];
		for (int i = 0; i < Lp; i++) {
			double dF_dPhi = 0;
			dF_dPhi = (phi_bar[i] + T*scikit.numerics.Math2.atanh(phi[i]));
			drift[i] = - dt*dF_dPhi;
		}

		double [] k2ft_drift = transform(drift);
		for (int i = 0; i < Lp; i++) {
			double kValue = (2*PI*i*R/L);
			k2ft_drift[2*i] *= (kValue*kValue);
			k2ft_drift[2*i+1] *= (kValue*kValue);
		}
		for (int i = 0; i < Lp; i++){
			phi_k[2*i] += k2ft_drift[2*i];
			phi_k[2*i+1] += k2ft_drift[2*i+1];
		}
		phi = backtransform (phi_k);
		
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
	private double [] transform(double [] src){

		for (int i = 0; i < Lp; i++) {
			fftScratch[2*i+0] = src[i];
			fftScratch[2*i+1] = 0;
		}
		fft.transform(fftScratch);
		return fftScratch;
		
	}


	private double [] backtransform(double [] src){

		fft.backtransform(src);
		double [] dest = new double [Lp];
		for (int i = 0; i < Lp; i++)
			dest[i] = src[2*i] / (Lp);
		return dest;
		
	}
	
}

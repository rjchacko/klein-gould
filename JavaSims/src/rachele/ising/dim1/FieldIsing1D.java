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
	public double t, noiseParam, noiseTerm, DENSITY;
	public double[] phi, F, phi_k, phi2, phi2_k;
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
		phi_bar = new double[Lp];
		del_phi = new double[Lp];
		phi2 = new double [Lp];
		phi2_k = new double[2*Lp];
		F = new double [Lp];
		phi_k = new double [2*Lp];
		
		fftScratch = new double[2*Lp];
		fft = new ComplexDoubleFFT_Mixed(Lp);
		
		for (int i = 0; i < Lp; i++){
			double newPhi = DENSITY + noiseParam*random.nextGaussian()*sqrt((1-DENSITY*DENSITY)/(dx*dx));
			phi[i] =  newPhi;
			phi2[i] = newPhi;
		}

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
			del_phi[i] = - dt*Lambda*(-phi_bar[i]-H+ T*atanh(phi[i]));
			del_phi[i] += sqrt(Lambda*dt*2*T/dx)*random.nextGaussian()*noiseParam;
		}
		for (int i = 0; i < Lp; i++)
			phi[i] += del_phi[i];	
		measureFreeEng();
		t += dt;
		phi_k = transform(phi);
	}

	public void simulateGlauber() {
		
		convolveWithRange(phi, phi_bar, R);	
		for (int i = 0; i < Lp; i++){
			double arg = (-phi_bar[i]/T + H/T);
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
			dF_dPhi = (-phi_bar[i] + T*scikit.numerics.Math2.atanh(phi[i]));
			drift[i] = - dt*dF_dPhi;
		}

		double [] k2ft_drift = transform(drift);

		for(int x = -Lp/2; x < Lp/2; x++){
			double kValue = (2*PI*x/Lp);
			int i = (x + Lp) % Lp;
			k2ft_drift[2*i] *= (kValue*kValue);
		}

		del_phi = backtransform(k2ft_drift);
		
		for (int i = 0; i < Lp; i++) {
			phi[i] += del_phi[i] + sqrt(dt*2*T/dx)*random.nextGaussian()*noiseParam;
		}	
		phi_k = transform(phi);
		t += dt;
	}

	public void simulateConservedFiniteDiff(){

		convolveWithRange(phi, phi_bar, R);
		double drift [] = new double [Lp];
		for (int i = 0; i < Lp; i++) {
			double dF_dPhi = 0;
			dF_dPhi = (-phi_bar[i] + T*scikit.numerics.Math2.atanh(phi[i]));
			drift[i] =  dt*dF_dPhi;
		}

		for (int i = 0; i < Lp; i++) {
			int irt = (i+1)/Lp;
			int ilf = (i-1+Lp)/Lp;
			del_phi[i] = (-2*drift[i] + drift[irt] + drift[ilf]);
		}
				
		for (int i = 0; i < Lp; i++) {
			phi[i] += del_phi[i] + sqrt(dt*2*T/dx)*random.nextGaussian()*noiseParam;
		}	
		phi_k = transform(phi);
		t += dt;
	}

	public void simulateConservedFiniteDiffMob(){

		convolveWithRange(phi, phi_bar, R);
		double drift [] = new double [Lp];
		for (int i = 0; i < Lp; i++) {
			double dF_dPhi = 0;
			dF_dPhi = (-phi_bar[i] + T*scikit.numerics.Math2.atanh(phi[i]));
			drift[i] =  dt*dF_dPhi;
		}

		double [] mgrad = new double [Lp];
		for (int i = 0; i < Lp; i++) {
			int irt = (i+1)%Lp;
			int ilf = (i-1+Lp)%Lp;
			int rrt = (i+2)%Lp;
			int llf = (i-2+Lp)%Lp;
			mgrad[i] = (-drift[rrt]+8*drift[irt]-8*drift[ilf]+drift[llf])/12;
		}
		
		
		for(int i = 0; i < Lp; i++){
			mgrad[i] *= (1-phi[i]*phi[i]);
		}
		
		for (int i = 0; i < Lp; i++) {
			int irt = (i+1)%Lp;
			int ilf = (i-1+Lp)%Lp;
			int rrt = (i+2)%Lp;
			int llf = (i-2+Lp)%Lp;
			del_phi[i] = (-mgrad[rrt]+8*mgrad[irt]-8*mgrad[ilf]+mgrad[llf])/12;
		}		
		
		for (int i = 0; i < Lp; i++) {	
			phi[i] += del_phi[i] + (1-phi[i]*phi[i])*sqrt(dt*2*T/dx)*random.nextGaussian()*noiseParam;
		}	
		phi_k = transform(phi);
		t += dt;
	}
	
	public void simulateConservedFspace(){

		convolveWithRange(phi, phi_bar, R);
		double drift [] = new double [Lp];
		for (int i = 0; i < Lp; i++) {
			double dF_dPhi = 0;
			dF_dPhi = (-phi_bar[i] + T*scikit.numerics.Math2.atanh(phi[i]));
			drift[i] = - dt*dF_dPhi;
		}

		double [] k2ft_drift = transform(drift);

		for(int x = -Lp/2; x < Lp/2; x++){
			double kValue = (2*PI*x/Lp);
			int i = (x + Lp) % Lp;
			k2ft_drift[2*i] *= (kValue*kValue);
		}


		for (int i = 0; i < Lp; i++) {
			phi_k[2*i] += k2ft_drift[2*i];
		}
		
		double [] scr = transform(phi);
		for (int i = 0; i < Lp; i++) {
			phi_k[2*i] = scr[2*i];
		}
		
		double [] scr2 = new double [2*Lp];
		for (int i = 0; i < 2*Lp; i++) {
			scr2[i] = phi_k[i];
		}

		double [] scr3 = backtransform(scr2);

		for (int i = 0; i < Lp; i++) {
			phi[i] = scr3[i]/2;
		}
		
		t += dt;
	}
	
	/**
	* Solves eqn of motion for system with conserved OP and 
	* mobility M=(1-phi(x)*phi(x))^2. 
	* 
	* This does not give the same effect as Glauber dynamics (ie- it still 
	* suffers from range exceptions with the arctanh.  Cannot use and additional
	* factor of M in real space, because this does not conserve mobility. For our 
	* purposes, it seems to work just the same as M=1 (program lasts the same amount
	* of time, so may as well just use M=1.)
	* 
	*  The eqn is:
	*  dphi/dt = grad dot mobility * grad(del F/ del phi)
	*  implimented as IFT(k*FT(mobility *IFT(k(FT(del F / del phi)))))
	*  or 
	*  k*f_(mob*{if_[k(f_dfdp)]})
	*/
	public void simulateConservedWithMobility(){
		convolveWithRange(phi, phi_bar, R);
		double drift [] = new double [Lp];
		for (int i = 0; i < Lp; i++) {
			double dF_dPhi = 0;
			dF_dPhi = (-phi_bar[i] + T*scikit.numerics.Math2.atanh(phi[i]));
			drift[i] = - dt*dF_dPhi;
		}

		double [] kft_drift = transform(drift);

		for(int x = -Lp/2; x < Lp/2; x++){
			double kValue = (2*PI*x/Lp);
			int i = (x + Lp) % Lp;
			kft_drift[2*i] *= (kValue);
		}
		fft.backtransform(kft_drift);
		fft.transform(kft_drift);
		for (int i = 0; i < Lp; i++) {
			kft_drift[2*i] *=  Math.pow(1-phi[i]*phi[i],2) / Lp;
			kft_drift[2*i+1] *= 1.0/Lp;
		}

		for(int x = -Lp/2; x < Lp/2; x++){
			double kValue = (2*PI*x/Lp);
			int i = (x + Lp) % Lp;
			kft_drift[2*i] *= (kValue);
		}		
		del_phi = backtransform(kft_drift);
		for (int i = 0; i < Lp; i++) {
			phi[i] += del_phi[i] + sqrt(dt*2*T/dx)*random.nextGaussian()*noiseParam;
		}	
		t += dt;
		System.out.println(del_phi[0]);
	}
	
	public void simulateConseveredSemiImp(){
		
		double [] atanh = new double [Lp];
		for (int i = 0; i < Lp; i++){
			atanh[i] = T*scikit.numerics.Math2.atanh(phi[i]);
		}
		double [] atanh_k = transform(atanh);
		for (int i = 0; i < Lp; i++) {
			fftScratch[2*i+0] = phi[i];
			fftScratch[2*i+1] = 0;
		}
		fft.transform(fftScratch);
		for (int i = 0; i < 2*Lp; i++) {
			phi_k[i] = fftScratch[i];
		}
		double [] del_phi_k = new double [2*Lp];
		for (int x = -Lp/2; x < Lp/2; x++) {
			double k = (2*PI*x/Lp);
			int i = (x + Lp) % Lp;
			double V = (k == 0 ? 1 : sin(k)/(k));
			double new_phi_k = (phi_k[2*i] + dt*k*k*atanh_k[2*i])/(1.0-k*k*J*V*dt);

			del_phi_k [2*i] = new_phi_k - phi_k[2*i];
//			System.out.println("del phi k " + k*k*atanh_k[2*i]);
//			new_phi_k = (phi_k[2*i+1] + dt*atanh_k[2*i+1])/(1.0-J*V*dt);
//			del_phi_k [2*i+1] = new_phi_k - phi_k[2*i+1];
		}
		del_phi = backtransform(del_phi_k);

		for (int i = 0; i < Lp; i++) {
			phi[i] += del_phi[i];// + sqrt(dt*2*T/dx)*random.nextGaussian()*noiseParam;
//			System.out.println("del phi " + del_phi_k[2*i]);
		}	
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
			del_phi[i] = - dt*Lambda*(-ampFactor*phi_bar[i]-H + T*atanh(phi[i]));
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
		for (int i = 0; i < Lp; i++){
			dest[i] = src[2*i] / (Lp);
		}
		return dest;	
	}
		
}

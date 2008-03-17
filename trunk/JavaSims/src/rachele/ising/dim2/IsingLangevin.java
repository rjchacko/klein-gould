package rachele.ising.dim2;

import static java.lang.Math.PI;
import static java.lang.Math.rint;
import static java.lang.Math.sin;
import static java.lang.Math.sqrt;
import static scikit.numerics.Math2.j1;
import static scikit.numerics.Math2.sqr;
import kip.util.Random;
import scikit.jobs.params.Parameters;
import scikit.numerics.fft.managed.ComplexDouble2DFFT;

public class IsingLangevin {
	public double L, R, T, dx, J, H;
	public int Lp, N;
	public double dt, t, maxTime;
	public double[] phi;
	double [] phi_bar, delPhi, Lambda, A;
	ComplexDouble2DFFT fft;
	double[] fftScratch;
	public static final double KR_SP = 5.13562230184068255630140;
	public static final double T_SP = 0.132279487396100031736846;	
	boolean noiselessDynamics = false;
	boolean circleInteraction = true;
	boolean magConservation = false;
	Random random = new Random();
	double DENSITY;
	
	public IsingLangevin(Parameters params) {
		random.setSeed(params.iget("Random seed", 0));
		J = params.fget("J");
		R = params.fget("R");
		L = R*params.fget("L/R");
		T = params.fget("T");
		H = params.fget("H");
		dx = R/params.fget("R/dx");
		dt = params.fget("dt");
		maxTime = params.fget("Max Time");
		DENSITY = params.fget("Magnetization");
		if(params.sget("Interaction") == "Circle") circleInteraction = true;
		else circleInteraction = false;
		if(params.sget("Noise") == "Off")noiselessDynamics = true;
		else noiselessDynamics = false;
		if(params.sget("Conserve M?") == "Yes") magConservation = true;
		else if(params.sget("Conserve M?") == "No") magConservation = false;
		Lp = Integer.highestOneBit((int)rint((L/dx)));
		dx = L / Lp;
		double RoverDx = R/dx;
		params.set("R/dx", RoverDx);
		N = Lp*Lp;
		t = 0;

		phi = new double[Lp*Lp];
		phi_bar = new double[Lp*Lp];
		delPhi = new double[Lp*Lp];
		Lambda = new double [Lp*Lp];
		fftScratch = new double[2*Lp*Lp];
		fft = new ComplexDouble2DFFT(Lp, Lp);
	}
	
	
	public void randomizeField(double m) {
		for (int i = 0; i < Lp*Lp; i++)
			phi[i] = m + random.nextGaussian()*sqrt((1-m*m)/(dx*dx));
	}
	
	public void readParams(Parameters params) {
		dt = params.fget("dt");
		H = params.fget("H");
		J = params.fget("J");
		R = params.fget("R");
		L = R*params.fget("L/R");
		dx = R/params.fget("R/dx");
		Lp = Integer.highestOneBit((int)rint((L/dx)));
		dx = L / Lp;
			
		params.set("R/dx", R/dx);
			
		if(params.sget("Interaction") == "Circle"){
			circleInteraction = true;
		}else{
			circleInteraction = false;
		}

		if(params.sget("Noise") == "Off"){
			noiselessDynamics = true;
		}else{
			noiselessDynamics = false;
		}

		if(params.sget("Conserve M?") == "Yes") magConservation = true;
		else if(params.sget("Conserve M?") == "No") magConservation = false;
	}
	
	void convolveWithRange(double[] src, double[] dest, double R) {
		double V;
		for (int i = 0; i < Lp*Lp; i++) {
			fftScratch[2*i] = src[i];
			fftScratch[2*i+1] = 0;
		}
		
		fft.transform(fftScratch);
		for (int y = -Lp/2; y < Lp/2; y++) {
			for (int x = -Lp/2; x < Lp/2; x++) {
				int i = Lp*((y+Lp)%Lp) + (x+Lp)%Lp;
				if(circleInteraction == true){
					double kR = (2*PI*sqrt(x*x+y*y)/L) * R;
					V = (kR == 0 ? 1 : circlePotential(kR));
				}else{
					double k_xR = (2*PI*x/L)*R;
					double k_yR =(2*PI*y/L)*R;
					V = (k_xR == 0 ? 1 : sin(k_xR)/k_xR);
					V *= (k_yR == 0 ? 1 : sin(k_yR)/k_yR);
				}
				fftScratch[2*i] *= V*J;
				fftScratch[2*i+1] *= V*J;
			}
		}
		fft.backtransform(fftScratch);
		
		for (int i = 0; i < Lp*Lp; i++) {
			dest[i] = fftScratch[2*i] / (Lp*Lp);
		}		
	}
	
	public void simulate() {
		convolveWithRange(phi, phi_bar, R);
		double meanLambda = 0;
		for (int i = 0; i < Lp*Lp; i++) {
			double dF_dPhi = 0;
			//dF_dPhi = -phi_bar[i]+T*(-log(1.0-phi[i])+log(1.0+phi[i]))/2.0 - H;
			dF_dPhi = -phi_bar[i]+T* scikit.numerics.Math2.atanh(phi[i])- H;
			Lambda[i] = 1; //shouldn't need lambda for linear theory
			delPhi[i] = - dt*Lambda[i]*dF_dPhi + sqrt(Lambda[i]*(dt*2*T)/dx)*noise();
			meanLambda += Lambda[i];
		}
		meanLambda /= Lp*Lp;
		double mu = (mean(delPhi)-(DENSITY-mean(phi)))/meanLambda;
		mu /= dt;
		if (magConservation == true){
			for (int i = 0; i < Lp*Lp; i++) 
				phi[i] += delPhi[i]-Lambda[i]*mu*dt;
		}
	}
	

	
	public void useNoiselessDynamics() {
		noiselessDynamics = true;
	}	
	
	public double phiVariance() {
		double var = 0;
		for (int i = 0; i < Lp*Lp; i++)
			var += sqr(phi[i]-DENSITY);
		return var / (Lp*Lp);
	}
	
	public double meanPhi(){
		double sum = 0;
		for (int i = 0; i < N; i ++)
			sum += phi[i];
		return sum/N;
	}
	
	public double squarePotential(double kRx, double kRy){
		return sin(kRx)*sin(kRy)/(kRx*kRy);
	}
	
	public double circlePotential(double kR){
		return 2*j1(kR)/kR;
	}
	
	double noise() {
		return noiselessDynamics ? 0 : random.nextGaussian();
	}
	
	public double mean(double[] a) {
		double sum = 0;
		for (int i = 0; i < Lp*Lp; i++)
				sum += a[i];
		return sum/(Lp*Lp); 
	}
//		
//	double meanSquared(double[] a) {
//		double sum = 0;
//		for (int i = 0; i < Lp*Lp; i++)
//				sum += a[i]*a[i];
//		return sum/(Lp*Lp);
//	}
	
	public double[] coarseGrained() {
		return phi;
	}
	
	public int numColumns() {
		return Lp;
	}
	
//	public double time() {
//		return t;
//	}
	
	public boolean circleInt(){
		return circleInteraction;
	}
	
//	public PointSet getHslice(){
//		int y = (int) (horizontalSlice * Lp);
//		double slice[] = new double[Lp];
//		for (int x = 0; x < Lp; x++) {
//			slice[x] = phi[Lp*y + x];
//		}
//		return new PointSet(0, dx, slice);
//	}
//
//	public PointSet getVslice(){
//		int x = (int) (verticalSlice * Lp);
//		double slice[] = new double[Lp];
//		for (int y = 0; y < Lp; y++) {
//			slice[y] = phi[Lp*y + x];
//		}
//		return new PointSet(0, dx, slice);
//	}
//
//	public PointSet get_delHslice(){
//		int y = (int) (horizontalSlice * Lp);
//		double slice[] = new double[Lp];
//		for (int x = 0; x < Lp; x++) {
//			slice[x] = delPhi[Lp*y + x];
//		}
//		return new PointSet(0, dx, slice);
//	}	
//
//	public PointSet get_delVslice(){
//		int x = (int) (verticalSlice * Lp);
//		double slice[] = new double[Lp];
//		for (int y = 0; y < Lp; y++) {
//			slice[y] = delPhi[Lp*y + x];
//		}
//		return new PointSet(0, dx, slice);
//	}	
}

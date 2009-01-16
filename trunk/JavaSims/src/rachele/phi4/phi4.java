package rachele.phi4;


import static java.lang.Math.rint;
import static java.lang.Math.sqrt;
import rachele.util.FourierTransformer;
import scikit.jobs.params.Parameters;
import scikit.numerics.fft.managed.ComplexDouble2DFFT;
import kip.util.Random;

public class phi4 {
	public double L, R, T, dx, H;
	public int Lp;
	
	public double J, dT;
	public int N;
	public double DENSITY;
	public double dt, t;
	public double[] phi, phiVector;
	public double freeEnergy;
	public double [] phi_bar, delPhi, Lambda, A;
	ComplexDouble2DFFT fft;	// Object to perform transforms
	FourierTransformer myfft;
	double[] fftScratch;
	public double [] mobility;  //spatially dependent mobility in general
	double noiseParameter;
	public boolean magConservation;
	Random random = new Random();
	
	public phi4(Parameters params) {
		random.setSeed(params.iget("Random seed", 1));
		R = params.fget("R");
		dx = R/params.fget("R/dx");
		L = R*params.fget("L/R");
		T = params.fget("T");
		dt = params.fget("dt");
		H = params.fget("H");
		DENSITY = params.fget("Magnetization");
		
		noiseParameter = params.fget("Noise");
		Lp = Integer.highestOneBit((int)rint((L/dx)));
		dx = L / Lp;
		System.out.println("dx = " + dx);
		double RoverDx = R/dx;
		params.set("R/dx", RoverDx);
		N = Lp*Lp;

		t = 0;

		phi = new double[Lp*Lp];
		phi_bar = new double[Lp*Lp];
		delPhi = new double[Lp*Lp];
		Lambda = new double [Lp*Lp];
		phiVector = new double[Lp*Lp];
		A = new double [Lp*Lp];
		mobility = new double [Lp*Lp];
		fftScratch = new double[2*Lp*Lp];
		fft = new ComplexDouble2DFFT(Lp, Lp);
		myfft = new FourierTransformer(Lp);
		randomizeField(DENSITY);
	}
	
	public void randomizeField(double m) {
		for (int i = 0; i < Lp*Lp; i++)
			
			//phi[i] = 0.00001*random.nextGaussian()/(dx);
			phi[i] = m + random.nextGaussian()*sqrt((1-m*m)/(dx*dx));
			//whether or not you start the system at the specified density doesn't matter too much because
			//the mag conservation algorithm will quickly take it to the correct density
	}
	
	public void readParams(Parameters params) {
		dt = params.fget("dt");
		R = params.fget("R");
		L = R*params.fget("L/R");
		dx = R/params.fget("R/dx");
		Lp = Integer.highestOneBit((int)rint((L/dx)));
		dx = L / Lp;
		H = params.fget("H");
		T = params.fget("T");
		params.set("R/dx", R/dx);
		params.set("Lp", Lp);
		noiseParameter = params.fget("Noise");

	}

	public void simulatePhi4(){
		//phi 4th theory from corberi paper
	
		for (int i = 0; i < Lp*Lp; i++){
			int x = i%Lp;
			int y = i/Lp;
			int rt = y*Lp + (x+1)%Lp;
			int lf = y*Lp + (x-1+Lp)%Lp;
			int up = Lp*((y + 1)%Lp) + x;
			int dn = Lp*((y - 1 + Lp)%Lp)+ x;
			
			double grad2 = -4*phi[i] + phi[rt] + phi[lf] + phi[up] + phi[dn];
			double g = 1.0;
			double gamma = 1.0;
			double r = -1.0;
//			double drift = dt*-gamma*(-R*grad2);
			double drift = dt*-gamma*(-R*grad2 + r*phi[i] +g*Math.pow(phi[i], 3));
//			double noiseTerm = noise()*sqrt(dt*2/(dx*dx));
			delPhi[i] =  drift;// + noiseTerm;
		}
		for (int i = 0; i < Lp*Lp; i++) 
			phi[i] += delPhi[i];			
		t += dt;		
	}
		
	public void restartClock(){t=0.0;}

		
	double noise() {
		//return noiselessDynamics ? 0 : random.nextGaussian();
		return noiseParameter*random.nextGaussian();
	}
	
	public double mean(double[] a) {
		double sum = 0;
		for (int i = 0; i < Lp*Lp; i++)
				sum += a[i];
		return sum/(Lp*Lp); 
	}
		
	double meanSquared(double[] a) {
		double sum = 0;
		for (int i = 0; i < Lp*Lp; i++)
				sum += a[i]*a[i];
		return sum/(Lp*Lp);
	}
	
	public double[] coarseGrained() {
		return phi;
	}
	
	public double time() {
		return t;
	}
		
	public void randomizeSeed(int newSeed){
		random.setSeed(newSeed);		
	}


}

package rachele.ising.dim2;

import static java.lang.Math.*;
import static scikit.numerics.Math2.*;
import java.io.DataInputStream;
import java.io.EOFException;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import kip.util.Random;
import rachele.util.FileUtil;
import rachele.util.FourierTransformer;
import scikit.dataset.Accumulator;
import scikit.dataset.PointSet;
import scikit.jobs.params.Parameters;
import scikit.numerics.fft.managed.ComplexDouble2DFFT;
import scikit.util.DoubleArray;

public class IsingField2D extends AbstractIsing2D{
	public double J, dT;
	public int N;
	public double DENSITY;
	public double dt, t;
	public double lastMu;
	public double[] phi, phiVector;
	public double freeEnergy;
	public int aveCount;
	public double [] phi_bar, delPhi, Lambda, A;
	ComplexDouble2DFFT fft;	// Object to perform transforms
	FourierTransformer myfft;
	double[] fftScratch;
	public double lambda;
	static double delta = 1e-10;
	public double switchD = 1;
	static double ZEPS = 0.0000000001;  //small number that protects against
	//trying to achieve fractional accuracy for a minimum that happens to be
	// exactly zero
	static double sigma = .2;
	public double [] mobility;  //spatially dependent mobility in general
	public double [] lineminDirection;
	public double [] xi;
	public double [] g;
	public double [] h;
	public double f_p;
	public double fret;
	public double driftContrib, absDriftContrib, absNoiseContrib;
	public double [] noiseTermC;
	public double [] driftTermC;
	
	Accumulator accClumpFreeEnergy;
	Accumulator accStripeFreeEnergy;
	Accumulator accEitherFreeEnergy;
	public Accumulator accFreeEnergy;
	double noiseParameter, stripeStrength;
	public boolean circleInteraction = false;
	public boolean magConservation;
	int slowPower = 0;
	String theory;
	Random random = new Random();
	
	public IsingField2D(Parameters params) {
		random.setSeed(params.iget("Random seed", 1));
		J = params.fget("J");
		R = params.fget("R");
		dx = R/params.fget("R/dx");
		L = R*params.fget("L/R");
		T = params.fget("T");
		dt = params.fget("dt");
		H = params.fget("H");
		//dT = params.fget("dT");
		//double dT = params.fget("dT");
		DENSITY = params.fget("Magnetization");
		accStripeFreeEnergy = new Accumulator();
		accClumpFreeEnergy = new Accumulator();
		accEitherFreeEnergy = new Accumulator();
		accFreeEnergy = new Accumulator(dt);
		
		if(params.sget("Interaction") == "Circle")
			circleInteraction = true;
		noiseParameter = params.fget("Noise");
		slowPower = 2;
		if(params.sget("Dynamics?") == "Langevin Conserve M") magConservation = true;
		else if(params.sget("Dynamics?") == "Langevin No M Conservation") magConservation = false;
		theory = params.sget("Approx");

		//theory = params.sget("Approx");
		Lp = Integer.highestOneBit((int)rint((L/dx)));
		dx = L / Lp;
		System.out.println("dx = " + dx);
		double RoverDx = R/dx;
		params.set("R/dx", RoverDx);
		//params.set("Lp", Lp);
		N = Lp*Lp;
		lineminDirection = new double [N];
		xi = new double [N];
		g = new double [N];
		h = new double [N];

		t = 0;

		phi = new double[Lp*Lp];
		phi_bar = new double[Lp*Lp];
		delPhi = new double[Lp*Lp];
		Lambda = new double [Lp*Lp];
		phiVector = new double[Lp*Lp];
		A = new double [Lp*Lp];
		mobility = new double [Lp*Lp];
		noiseTermC = new double [Lp*Lp];
		driftTermC = new double [Lp*Lp];
		//HH = new double[Lp*Lp];
		//double backgroundH = params.fget("H");
		//double stripeH = params.fget("H Stripe");
		//setExternalField(backgroundH, stripeH);
		
		fftScratch = new double[2*Lp*Lp];
		fft = new ComplexDouble2DFFT(Lp, Lp);
		myfft = new FourierTransformer(Lp);
		
	
		//String init = "Random Gaussian";
		String init = params.sget("Init Conditions");
		if(init == "Random Gaussian"){
			randomizeField(DENSITY);
			for (int i = 0; i < Lp*Lp; i++)
				mobility[i] = 1.0;
//				mobility[i] = (1.0 - DENSITY*DENSITY)/T;
				System.out.println("Random Gaussian");
		}else if(init == "Read 1D Soln"){
			System.out.println("Read in 1D solution");
			String fileName = params.sget("1D Input File");
			double[] phiSlice = new double [Lp];  
			phiSlice = getSymmetricSlice(fileName);
			set1DConfig(phiSlice);
			for (int i = 0; i < Lp*Lp; i++)
				mobility[i] = 1;
//				mobility[i] = (1 - phiSlice[i%Lp]*phiSlice[i%Lp]) / T;
		}else if (init == "Constant"){
			System.out.println("Constant");
			for (int i = 0; i < Lp*Lp; i ++)
				phi[i] = DENSITY;
			for (int i = 0; i < Lp*Lp; i++)
				mobility[i] = (1.0 - DENSITY*DENSITY)/T;
		}else if(init == "Read From File"){
			readInitialConfiguration(params.sget("2D Input Dir") + File.separator + "inputConfig");
			System.out.println("Read From File");
			for (int i = 0; i < Lp*Lp; i++)
				mobility[i] = (1.0 - DENSITY*DENSITY)/T;			
		}else if(init == "Artificial Stripe 3"){
			System.out.println("Artificial Stripe 3");
			int stripeLoc = (int)(Lp/3.0);
			randomizeField(DENSITY);
			double m = DENSITY + (1.0-DENSITY)*stripeStrength;
			for (int i = 0; i < Lp*Lp; i ++){
				int x = i%Lp;
				if (x == stripeLoc)
					phi[i] = m + random.nextGaussian()*sqrt((1-m*m)/(dx*dx));
				if (x == stripeLoc*2)
					phi[i] = m + random.nextGaussian()*sqrt((1-m*m)/(dx*dx));
				if (x == stripeLoc*3)
					phi[i] = m + random.nextGaussian()*sqrt((1-m*m)/(dx*dx));
				if (phi[i] > 1.0){
					System.out.println(i + " " + phi[i]);
				}
			}
			for (int i = 0; i < Lp*Lp; i++)   //This is just filler.  Constant mobility in this situation doesn't mean much.
				mobility[i] = (1.0 - DENSITY*DENSITY)/T;
		}else if(init == "Artificial Stripe 2"){
			System.out.println("Artificial Stripe 2");
			randomizeField(DENSITY);
			double m = DENSITY - (1.0-DENSITY)*stripeStrength;
			for (int i = 0; i < Lp; i ++)
				phi[i] = m + random.nextGaussian()*sqrt((1-m*m)/(dx*dx));
			for (int i = Lp*Lp/2; i < Lp + Lp*Lp/2; i ++)
				phi[i] = m + random.nextGaussian()*sqrt((1-m*m)/(dx*dx));
			for (int i = 0; i < Lp*Lp; i++)   //This is just filler.  Constant mobility in this situation doesn't mean much.
				mobility[i] = (1.0 - DENSITY*DENSITY)/T;		
		}
	}

	public void set1DConfig(double [] config){
		for(int i = 0; i < Lp*Lp; i++)
//			phi[i]=config[i%Lp] + (random.nextGaussian())*noiseParameter*sqrt(1-pow(config[i%Lp],2))/(dx*dx);
		phi[i]=config[i%Lp] + (random.nextGaussian())*sqrt(1-pow(config[i%Lp],2))/(dx*dx);
	}

	public double [] getSymmetricSlice(String file){
		double [] slice = new double [Lp];
		double [] temp= FileUtil.readConfigFromFile(file, Lp);
		double minPhi0Value = 1.0;
		int minPhi0Location = -1;
		for (int i = 0; i < Lp; i++){
			if (temp[i] < minPhi0Value){
				minPhi0Location = i;
				minPhi0Value = temp[i];
			}
		}	
		for (int i = 0; i < Lp; i++){
			slice[i] = temp[(minPhi0Location+i)%Lp];
		}
		return slice;
	}


	public void readInitialConfiguration(String fileName){
		try{
			File myFile = new File(fileName);
			DataInputStream dis = new DataInputStream(new FileInputStream(myFile));
			int spaceIndex;
			double phiValue;
			try{
				while(true){
					spaceIndex =dis.readInt();
					dis.readChar();       // throws out the tab
					phiValue = dis.readDouble();
					dis.readChar();
					phi[spaceIndex] = phiValue;
				}
			} catch (EOFException e) {
			}

		} catch (FileNotFoundException e) {
			System.err.println("FileStreamsTest: " + e);
		} catch (Exception ex) {
			ex.printStackTrace();
		}
	}
	
	public void randomizeField(double m) {
		for (int i = 0; i < Lp*Lp; i++)
			
			//phi[i] = 0.00001*random.nextGaussian()/(dx);
			phi[i] = m + random.nextGaussian()*sqrt((1-m*m)/(dx*dx));
			//whether or not you start the system at the specified density doesn't matter too much because
			//the mag conservation algorithm will quickly take it to the correct density
	}
	
	public void readParams(Parameters params) {
		//if (params.sget("Plot FEvT") == "Off") T = params.fget("T");
		dt = params.fget("dt");
		J = params.fget("J");
		R = params.fget("R");
		L = R*params.fget("L/R");
		dx = R/params.fget("R/dx");
		Lp = Integer.highestOneBit((int)rint((L/dx)));
		dx = L / Lp;
		H = params.fget("H");
		T = params.fget("T");
		params.set("R/dx", R/dx);
		params.set("Lp", Lp);
	
		if(params.sget("Interaction") == "Circle") circleInteraction = true;
		else circleInteraction = false;
		
		noiseParameter = params.fget("Noise");
		if(params.sget("Dynamics?") == "Langevin Conserve M") magConservation = true;
		else if(params.sget("Dynamics?") == "Langevin No M Conservation") magConservation = false;
		theory = params.sget("Approx");
	}
	
	public void initializeFieldWithSeed() {
		for (int i = 0; i < Lp*Lp; i++) {
			double x = dx*(i%Lp - Lp/2);
			double y = dx*(i/Lp - Lp/2);
			double r = sqrt(x*x+y*y);
			double mag = 0.8 / (1+sqr(r/R));
			
			double kR = KRcircle; // it's fun to try different values
			double x1 = x*cos(1*PI/6) + y*sin(1*PI/6);
			double x2 = x*cos(3*PI/6) + y*sin(3*PI/6);
			double x3 = x*cos(5*PI/6) + y*sin(5*PI/6);
			phi[i] = DENSITY*(1+mag*(cos(x1*kR/R) + cos(x2*kR/R) + cos(x3*kR/R)));			
		}
	}
	
	public void convolveWithRange(double[] src, double[] dest, double R) {
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
					V = (kR == 0 ? 1 : 2*j1(kR)/kR);
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

	public void simulateSimple(){
		//looks like its working
		convolveWithRange(phi, phi_bar, R);	
		for (int i = 0; i < Lp*Lp; i++){
			delPhi[i] = -dt*mobility[i]*(-phi_bar[i]+T* scikit.numerics.Math2.atanh(phi[i])- H)+ noise()*sqrt(dt*2/(dx*dx));
		}
		for (int i = 0; i < Lp*Lp; i++) 
			phi[i] += delPhi[i];			
		t += dt;		
	}
	
	public void simulateGlauber(){
		convolveWithRange(phi, phi_bar, R);	
		for (int i = 0; i < Lp*Lp; i++){
			double arg = (phi_bar[i]/T + H/T);
			double tanh = Math.tanh(arg);
			double driftTerm = dt*(tanh-phi[i]);
			double noisePre = sqrt(2-pow(Math.tanh(arg),2)-pow(phi[i],2));
			double noiseTerm = noisePre*noise()*sqrt(dt*2/(dx*dx));
			delPhi[i] = driftTerm + noiseTerm;
		}

		for (int i = 0; i < Lp*Lp; i++) 
			phi[i] += delPhi[i];			
		t += dt;	
	}

	public void glauberCalcContrib() {
		convolveWithRange(phi, phi_bar, R);
		driftContrib = 0.0;
		absDriftContrib = 0.0;
		absNoiseContrib = 0.0;

		for (int i = 0; i < Lp*Lp; i++){
			double arg = (phi_bar[i] + H)/T;
			double tanh = Math.tanh(arg);
			driftTermC[i] = dt*(tanh-phi[i]);
			double noisePre = sqrt(2-pow(Math.tanh(arg),2)-pow(phi[i],2));
			noiseTermC[i] = noisePre*noise()*sqrt(dt*2/(dx*dx));
			delPhi[i] = driftTermC[i] + noiseTermC[i];
			//driftContrib += abs(tanh)/(abs(tanh)+abs(phi[i]));
			driftContrib += abs(driftTermC[i])/(abs(driftTermC[i])+abs(noiseTermC[i]));
			absDriftContrib += abs(driftTermC[i]);
			absNoiseContrib += abs(noiseTermC[i]);
		}
		for (int i = 0; i < Lp*Lp; i++) 
			phi[i] += delPhi[i];			
		driftContrib /= (Lp*Lp);
		absNoiseContrib /= (Lp*Lp);
		absDriftContrib /= (Lp*Lp);
		t += dt;	
	}
	
	
	
	public boolean simulateUnstable(){
		boolean abortMode = true;
		boolean abort = false; //abort simulation with range exception?
		convolveWithRange(phi, phi_bar, R);		
		double dt0 = dt;
		for (int i = 0; i < Lp*Lp; i++){
			delPhi[i] = -dt0*mobility[i]*(-phi_bar[i]+T* scikit.numerics.Math2.atanh(phi[i])- H)+ sqrt(dt0*2/(dx*dx))*noise();
		}
		double [] testPhi = new double [Lp*Lp];
		if(abortMode){
			for (int i = 0; i < Lp*Lp; i++){
				testPhi[i] = phi[i]+ delPhi[i];
				if(testPhi[i] >= 1.0){ 
					abort = true;
				}else if(testPhi[i] <= -1.0){
					abort = true;
				}
			}
			if(abort){
				System.out.println("abort at " + time());
			}else
				for (int i = 0; i < Lp*Lp; i++) phi[i] += dt*delPhi[i]/dt0;			
			t += dt;
		}else{//older code:  keeps the program going but kind of screws it up
			int ct = 0;
			for (int i = 0; i < Lp*Lp; i++){
				testPhi[i] = phi[i]+ delPhi[i];
				if(testPhi[i] > 1.0){
					ct += 1;
					double dtTest = dt0*((1-phi[i])/2.0)/testPhi[i];
					if(dtTest < dt) dt = dtTest;
				}else if(testPhi[i] < -1.0){
					ct += 1;
					double dtTest = dt0*((-1-phi[i])/2.0)/testPhi[i];
					if(dtTest < dt) dt = dtTest;
				}
			}
			if (ct != 0) System.out.println("ct = " + ct);
		
			for (int i = 0; i < Lp*Lp; i++) 
				phi[i] += dt*delPhi[i]/dt0;			
			t += dt;
			//find max, min and average delphis	
			double meanDelPhi = mean(delPhi);
			double maxDelPhi = DoubleArray.max(delPhi);
			double minDelPhi = DoubleArray.min(delPhi);
			System.out.println("mean = " + meanDelPhi + " min = " + minDelPhi + " max = " + maxDelPhi);

		}
		return abort;
	}
	
	public void restartClock(){t=0.0;}

	/**
	* This purpose of this method is to simulate with modified
	* dynamics and calculate the average percent contribution 
	* from the drift term (dF/dphi term) and the noise term.
	* Then should be able to test with zero init noise which
	* term "drives" the dynamics and where does it switch from one 
	* term to another. The main reason for doing this is to 
	* investigate the disorder->clump->stripe transition.
	*/
	public void simModCalcContrib() {
		convolveWithRange(phi, phi_bar, R);
		driftContrib = 0.0;
		for (int i = 0; i < Lp*Lp; i++) {
			double dF_dPhi = 0;
			dF_dPhi = mobility[i]*(-phi_bar[i] - H + T*scikit.numerics.Math2.atanh(phi[i]));
			Lambda[i] = sqr(1 - phi[i]*phi[i]);		
			double driftTerm = - dt*Lambda[i]*dF_dPhi;
			double noiseTerm = sqrt(Lambda[i]*(dt*2*T*mobility[i])/dx)*noise();
			delPhi[i] = driftTerm + noiseTerm;
			//System.out.println(noiseTerm + " " + driftTerm + " " + delPhi[i]);
			driftContrib += abs(driftTerm)/(abs(driftTerm)+abs(noiseTerm));
			phiVector[i] = delPhi[i];
		}
		for (int i = 0; i < Lp*Lp; i++) 
			phi[i] += delPhi[i];	
		driftContrib /= (Lp*Lp);
		t += dt;
	}

	public void simCalcContrib() {
		convolveWithRange(phi, phi_bar, R);
		driftContrib = 0.0;
		for (int i = 0; i < Lp*Lp; i++) {
			double dF_dPhi = 0;
			dF_dPhi = mobility[i]*(-phi_bar[i] - H + T*scikit.numerics.Math2.atanh(phi[i]));
			double driftTerm = - dt*dF_dPhi;
			double noiseTerm = sqrt((dt*2*T*mobility[i])/dx)*noise();
			delPhi[i] = driftTerm + noiseTerm;
			//System.out.println(noiseTerm + " " + driftTerm + " " + delPhi[i]);
			driftContrib += abs(driftTerm/dt)/(abs(driftTerm/dt)+abs(noiseTerm/dt));
			phiVector[i] = delPhi[i];
		}
		for (int i = 0; i < Lp*Lp; i++) 
			phi[i] += delPhi[i];	
		driftContrib /= (Lp*Lp);
		t += dt;
	}
	
	public double [] simModCalcContribK(int kyValue) {
		convolveWithRange(phi, phi_bar, R);
		double [] driftProp = new double [Lp];
		for (int i = 0; i < Lp*Lp; i++) {
			double dF_dPhi = 0;
			dF_dPhi = mobility[i]*(-phi_bar[i] - H + T*scikit.numerics.Math2.atanh(phi[i]));
			Lambda[i] = sqr(1 - phi[i]*phi[i]);		
			double driftTerm = - dt*Lambda[i]*dF_dPhi;
			double noiseTerm = sqrt(Lambda[i]*(dt*2*T*mobility[i])/dx)*noise();
			delPhi[i] = driftTerm + noiseTerm;
			//System.out.println(noiseTerm + " " + driftTerm + " " + delPhi[i]);
			driftContrib += abs(driftTerm)/(abs(driftTerm)+abs(noiseTerm));
			phiVector[i] = delPhi[i];
			if (i/Lp == kyValue) driftProp[i%Lp] = driftContrib;
		}
		for (int i = 0; i < Lp*Lp; i++) 
			phi[i] += delPhi[i];	
		driftContrib /= (Lp*Lp);
		t += dt;
		return driftProp;
		
	}

	
	public void simulateConserved(){
		convolveWithRange(phi, phi_bar, R);
		double drift [] = new double [Lp*Lp];
		for (int i = 0; i < Lp*Lp; i++) {
			double dF_dPhi = 0;
			dF_dPhi = (-phi_bar[i] + T*scikit.numerics.Math2.atanh(phi[i]));
			drift[i] = - dt*dF_dPhi;
		}

		double [] kft_drift = transform(drift);
		for (int i = 0; i < Lp*Lp; i++) {
			int ky = i / Lp;
			int kx = i % Lp;
			double kValue = (2*PI*sqrt((double)(ky*ky+kx*kx))*R/L);
			kft_drift[2*i] *= (kValue);
			kft_drift[2*i+1] *= (kValue);
		}
		double [] m_kdrift = backtransform(kft_drift);
		for (int i = 0; i < Lp*Lp; i++) 
			m_kdrift[i] *= Math.pow(1-phi[i]*phi[i],2);
		
		double [] k_mkdrift = transform(m_kdrift);
		for (int i = 0; i < Lp*Lp; i++) {
			int ky = i / Lp;
			int kx = i % Lp;
			double kValue = (2*PI*sqrt((double)(ky*ky+kx*kx))*R/L);
			k_mkdrift[2*i] *= (kValue);
			k_mkdrift[2*i+1] *= (kValue);
		}		
		
		delPhi = backtransform(k_mkdrift);
		
		for (int i = 0; i < Lp*Lp; i++) {
			double noiseTerm = sqrt(dt*2/(dx*dx))*noise();
			phi[i] += delPhi[i]+noise()*noiseTerm;
			driftTermC[i] = delPhi[i];
			driftContrib += abs(driftTermC[i])/(abs(driftTermC[i])+abs(noiseTerm));
			absDriftContrib += abs(driftTermC[i]);
			absNoiseContrib += abs(noiseTerm);
			
		}	
		driftContrib /= (Lp*Lp);
		absNoiseContrib /= (Lp*Lp);
		absDriftContrib /= (Lp*Lp);
		t += dt;
//		System.out.println("conserved");

	}
	
	public void simulate() {
		freeEnergy = 0;  //free energy is calculated for previous time step
		double potAccum = 0;
		double entAccum = 0;
		double del_phiSquared = 0;
		//theory = "Exact";//use exact for halfStep rules- halfStep does not allow phi to change more than 1/2 way to singularity
		
		convolveWithRange(phi, phi_bar, R);
		
		double meanLambda = 0;
		
		for (int i = 0; i < Lp*Lp; i++) {
			double dF_dPhi = 0, entropy = 0;
			if (theory == "Phi4" || theory == "Phi4HalfStep"){
				dF_dPhi = mobility[i]*(-phi_bar[i]+T*(phi[i]+pow(phi[i],3)) - H);	
				Lambda[i] = 1;
			}else{
				dF_dPhi = mobility[i]*(-phi_bar[i] - H + T*scikit.numerics.Math2.atanh(phi[i]));
				//dF_dPhi = -phi_bar[i]+T*(-log(1.0-phi[i])+log(1.0+phi[i]))/2.0 - H;
				if(theory == "HalfStep"){
					Lambda[i] = 1;
					entropy = -((1.0 + phi[i])*log(1.0 + phi[i]) +(1.0 - phi[i])*log(1.0 - phi[i]))/2.0;
				}else{
					Lambda[i] = sqr(1 - phi[i]*phi[i]);				
					entropy = -((1.0 + phi[i])*log(1.0 + phi[i]) +(1.0 - phi[i])*log(1.0 - phi[i]))/2.0;
				}
			}
			delPhi[i] = - dt*Lambda[i]*dF_dPhi + sqrt(Lambda[i]*(dt*2*T*mobility[i])/dx)*noise();
			phiVector[i] = delPhi[i];
			meanLambda += Lambda[i];
			double potential = -(phi[i]*phi_bar[i])/2.0;
			potAccum += potential;
			entAccum -= T*entropy;
			freeEnergy += potential - T*entropy - H*phi[i];

		}		
		meanLambda /= Lp*Lp;
		double mu = (mean(delPhi)-(DENSITY-mean(phi)))/meanLambda;
		mu /= dt;
		if (magConservation == false)
			mu = 0;
		for (int i = 0; i < Lp*Lp; i++) {
			freeEnergy +=  -mu*phi[i];
			delPhi[i] -= Lambda[i]*mu*dt;
			double testPhi = phi[i]+ delPhi[i];
			if(theory == "HalfStep" || theory == "Phi4HalfStep"){
				if(abs(testPhi) >= 1){
					if(testPhi>=1)
						testPhi=phi[i]+(1-phi[i])/2.0;
					else if(testPhi<=-1)
						testPhi=phi[i]+(-1-phi[i])/2.0;
				}
			}	
			//phi[i] += delPhi[i]-Lambda[i]*mu*dt;
			//phi[i] += delPhi[i];
			phi[i]=testPhi;
			del_phiSquared += phi[i]*phi[i];
		}
		lastMu = mu;
		freeEnergy /= (Lp*Lp) ;
		//accFreeEnergy.accum(t, freeEnergy);
		if(theory == "TimeAdjust") dt = mean(Lambda);
		t += dt;
	}
	
	public void changeT(){
		T += dT;
		System.out.println("T = " + T + " " + dT);
	}
	
	public void accumClumpFreeEnergy(){
		accClumpFreeEnergy.accum(T, freeEnergy);
	}
	
	public void accumStripeFreeEnergy(){
		accStripeFreeEnergy.accum(T, freeEnergy);
	}

	public void accumEitherFreeEnergy(){
		accEitherFreeEnergy.accum(T, freeEnergy);
	}
	
	public double phiVariance() {
		double var = 0;
		for (int i = 0; i < Lp*Lp; i++)
			var += sqr(phi[i]-DENSITY);
		return var / (Lp*Lp);
	}
	
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
	
	public int numColumns() {
		return Lp;
	}
	
	public double time() {
		return t;
	}
	
	public boolean circleInt(){
		return circleInteraction;
	}
	
	public PointSet getHslice(double horizontalSlice){
		int y = (int) (horizontalSlice * Lp);
		double slice[] = new double[Lp];
		for (int x = 0; x < Lp; x++) {
			slice[x] = phi[Lp*y + x];
		}
		return new PointSet(0, 1.0, slice);
	}

	public PointSet getVslice(double verticalSlice){
		int x = (int) (verticalSlice * Lp);
		double slice[] = new double[Lp];
		for (int y = 0; y < Lp; y++) {
			slice[y] = phi[Lp*y + x];
		}
		return new PointSet(0, 1.0, slice);
	}

	public PointSet get_delHslice(double horizontalSlice){
		int y = (int) (horizontalSlice * Lp);
		double slice[] = new double[Lp];
		for (int x = 0; x < Lp; x++) {
			slice[x] = delPhi[Lp*y + x];
		}
		return new PointSet(0, dx, slice);
	}	

	public PointSet get_delVslice(double verticalSlice){
		int x = (int) (verticalSlice * Lp);
		double slice[] = new double[Lp];
		for (int y = 0; y < Lp; y++) {
			slice[y] = delPhi[Lp*y + x];
		}
		return new PointSet(0, dx, slice);
	}	
	
	public Accumulator getFreeEnergyAcc() {
		return accFreeEnergy;
	}
	
	public Accumulator getStripeFreeEnergyAcc() {
		return accStripeFreeEnergy;
	}
	
	public Accumulator getClumpFreeEnergyAcc() {
		return accClumpFreeEnergy;
	}
	
	public Accumulator getEitherFreeEnergyAcc() {
		return accEitherFreeEnergy;
	}
	
	public void initializeConjGrad(){
    	// Initializations:
    	// evaluate function and derivative at given point
    	// set all vectors (xi, g, h) equal to the direction
    	// of steepest Descent at this point
		f_p = isingFreeEnergyCalc(phi);
		xi = steepestAscentCalc(phi);
    	for(int j = 0; j < N; j ++){
    		g[j] = -xi[j];
    		h[j] = g[j];
    		xi[j] = h[j];
    	}		
	}
	
	public double isingFreeEnergyCalc(double [] config){
		convolveWithRange(config, phi_bar, R);
		freeEnergy = 0;
		for (int i = 0; i < Lp*Lp; i++) {
			if(Math.abs(phi[i]) >= 1)
				System.out.println("boundary at point " + i);
			double entropy = -((1.0 + config[i])*log(1.0 + config[i]) +(1.0 - config[i])*log(1.0 - config[i]))/2.0;
			double potential = -(config[i]*phi_bar[i])/2.0;
			freeEnergy += potential - T*entropy - H*config[i];
		}
		freeEnergy /= (Lp*Lp);
		if (Double.isNaN(freeEnergy))
			return Double.POSITIVE_INFINITY;
		return freeEnergy;
	}

	public double [] steepestAscentCalc(double [] config){
		double steepestAscentDir [] = new double [N];
		convolveWithRange(config, phi_bar, R);
		for (int i = 0; i < Lp*Lp; i++) {
			steepestAscentDir[i] = (-phi_bar[i] +T* scikit.numerics.Math2.atanh(config[i])- H);
			steepestAscentDir[i] *= (pow(1 - phi[i]*phi[i], slowPower));
		}
		return steepestAscentDir;		
	}
	
	public double isingFreeEnergyCalcA (double [] config){
		convolveWithRange(config, phi_bar, R);
		freeEnergy = 0;
		for (int i = 0; i < Lp*Lp; i++) {
			if(Math.abs(phi[i]) >= 1)
				System.out.println("boundary at point " + i);
			double entropy = -((1.0 + config[i])*log(1.0 + config[i]) +(1.0 - config[i])*log(1.0 - config[i]))/2.0;
			double potential = -(config[i]*phi_bar[i])/2.0;
			//double boundaryTerm = exp(-sqr(1-config[i])/(2.0*sigma*sigma))/sqr(1-config[i]) + exp(-sqr(1+config[i])/(2.0*sigma*sigma))/sqr(1+config[i]);
			//System.out.println("Boundary term = " + boundaryTerm);
			freeEnergy += potential - T*entropy*100 - H*config[i];
		}
		freeEnergy /= (Lp*Lp);
		if (Double.isNaN(freeEnergy))
			return Double.POSITIVE_INFINITY;
		return freeEnergy;
	}

	public double [] steepestAscentCalcA(double [] config){
		double steepestAscentDir [] = new double [N];
		convolveWithRange(config, phi_bar, R);
		for (int i = 0; i < Lp*Lp; i++) {
			double boundaryTerm =  exp(-sqr(1-config[i])/(2.0*sigma*sigma))*(1/sqr(sigma) + 2.0 / (sqr(1-config[i])))/(1-config[i]);
			//boundaryTerm -=        exp(-sqr(1+config[i])/(2.0*sigma*sigma))*(1/sqr(sigma) + 2.0 / (sqr(1+config[i])))/(1+config[i]);
			//System.out.println("Boundary term = " + boundaryTerm);
			steepestAscentDir[i] = (-phi_bar[i] +T*100* scikit.numerics.Math2.atanh(config[i])- H + boundaryTerm);
			//steepestAscentDir[i] = (-phi_bar[i] +T* configA[i]- H + boundaryTerm);///(1-sqr(Math.tanh(configA[i])));//*(sqr(1 - phi[i]*phi[i]));
		}
		return steepestAscentDir;		
	}
	
	public void randomizeSeed(int newSeed){
		random.setSeed(newSeed);		
	}

	public void initFile(Parameters params, String fileName, String message1, String message2){
		FileUtil.printlnToFile(fileName, message1);
		FileUtil.printlnToFile(fileName, message2);
		FileUtil.printlnToFile(fileName, "# Parameters follow:");
		FileUtil.printlnToFile(fileName, "# Interaction = ", params.sget("Interaction"));
		FileUtil.printlnToFile(fileName, "# Dynamics = ", params.sget("Dynamics?"));
		FileUtil.printlnToFile(fileName, "# init = ", params.sget("Init Conditions"));
		FileUtil.printlnToFile(fileName, "# Approx = ", params.sget("Approx"));
		FileUtil.printlnToFile(fileName, "# Noise = ", params.sget("Noise"));
		FileUtil.printlnToFile(fileName, "# Random Seed = ", params.sget("Random seed"));
		FileUtil.printlnToFile(fileName, "# L/R = ", params.sget("L/R"));
		FileUtil.printlnToFile(fileName, "# R = ", params.sget("R"));
		FileUtil.printlnToFile(fileName, "# R/dx  = ", params.sget("R/dx"));
		FileUtil.printlnToFile(fileName, "# init mag = ", params.sget("Magnetization"));
		FileUtil.printlnToFile(fileName, "# temperature ", params.sget("T"));
		FileUtil.printlnToFile(fileName, "# Random Seed = ", params.sget("Random seed"));
		FileUtil.printlnToFile(fileName, "# J = ", params.sget("J"));
		FileUtil.printlnToFile(fileName, "# h = ", params.sget("H"));
		FileUtil.printlnToFile(fileName, "# dt = ", params.sget("dt"));
		FileUtil.printlnToFile(fileName, "# max time = ", params.sget("Max Time"));
		FileUtil.printlnToFile(fileName, "# reps = ", params.sget("Reps"));
	}
	private double [] transform(double [] src){

		for (int i = 0; i < Lp*Lp; i++) {
			fftScratch[2*i+0] = src[i];
			fftScratch[2*i+1] = 0;
		}
		fft.transform(fftScratch);
		return fftScratch;
	}


	private double [] backtransform(double [] src){
		fft.backtransform(src);
		double [] dest = new double [Lp*Lp];
		for (int i = 0; i < Lp*Lp; i++)
			dest[i] = src[2*i] / (Lp*Lp);
		return dest;
	}
//	private void setExternalField(double backgroundH, double stripeH){
//		for (int i = Lp; i < Lp*Lp; i++)
//			HH[i] = backgroundH;
//		for (int i = 0; i < Lp; i++)
//			HH[i] = stripeH;
//		for (int i = Lp*Lp/2; i < Lp + Lp*Lp/2; i++)
//			HH[i] = stripeH;	
//	}
}

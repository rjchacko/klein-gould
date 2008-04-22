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
import scikit.dataset.Accumulator;
import scikit.dataset.PointSet;
import scikit.jobs.params.Parameters;
import scikit.numerics.fft.managed.ComplexDouble2DFFT;

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
	public double horizontalSlice;
	public double verticalSlice;
	ComplexDouble2DFFT fft;	// Object to perform transforms
	double[] fftScratch;
	
	public double lambda;
	static double delta = 1e-10;
	public double switchD = 1;
	static double ZEPS = 0.0000000001;  //small number that protects against
	//trying to achieve fractional accuracy for a minimum that happens to be
	// exactly zero
	static double sigma = .2;

	public double [] lineminDirection;
	public double [] xi;
	public double [] g;
	public double [] h;
	
	public double f_p;
	public double fret;
	
	Accumulator accClumpFreeEnergy;
	Accumulator accStripeFreeEnergy;
	Accumulator accEitherFreeEnergy;
	public Accumulator accFreeEnergy;
			
	//boolean noiselessDynamics = false;
	double noiseParameter, stripeStrength;
	boolean circleInteraction = false;
	boolean magConservation = false;
	int slowPower = 0;
	String theory;
	
	Random random = new Random();
	

	
	public IsingField2D(Parameters params) {
		random.setSeed(params.iget("Random seed", 0));
		J = params.fget("J");
		R = params.fget("R");
		L = R*params.fget("L/R");
		T = params.fget("T");
		dx = R/params.fget("R/dx");
		dt = params.fget("dt");
		H = params.fget("H");
		//dT = params.fget("dT");
		//double dT = params.fget("dT");
		DENSITY = params.fget("Magnetization");
		accStripeFreeEnergy = new Accumulator(0.0001);
		accClumpFreeEnergy = new Accumulator(0.0001);
		accEitherFreeEnergy = new Accumulator(0.0001);
		accFreeEnergy = new Accumulator(dt);
		
		horizontalSlice = params.fget("Horizontal Slice");
		verticalSlice = params.fget("Vertical Slice");
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
		//HH = new double[Lp*Lp];
		//double backgroundH = params.fget("H");
		//double stripeH = params.fget("H Stripe");
		//setExternalField(backgroundH, stripeH);
		
		fftScratch = new double[2*Lp*Lp];
		fft = new ComplexDouble2DFFT(Lp, Lp);
		
		//String init = "Random Gaussian";
		String init = params.sget("Init Conditions");
		if(init == "Random Gaussian"){
			randomizeField(DENSITY);
			System.out.println("Random Gaussian");
		}else if(init == "Read 1D Soln"){
			System.out.println("Read in 1D solution");
			String fileName = "../../../research/javaData/configs1d/config";
			double[] phiSlice = new double [Lp];
			phiSlice = FileUtil.readConfigFromFile(fileName, Lp);
			for(int i = 0; i < Lp*Lp; i++)
				phi[i]=phiSlice[i%Lp] + (random.nextGaussian())/1000000;
		}else if (init == "Constant"){
			System.out.println("Constant");
			for (int i = 0; i < Lp*Lp; i ++)
				phi[i] = DENSITY;
		}else if(init == "Read From File"){
			readInitialConfiguration();
			System.out.println("Read From File");
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
		}else if(init == "Artificial Stripe 2"){
			System.out.println("Artificial Stripe 2");
			randomizeField(DENSITY);
			double m = DENSITY - (1.0-DENSITY)*stripeStrength;
			for (int i = 0; i < Lp; i ++)
				phi[i] = m + random.nextGaussian()*sqrt((1-m*m)/(dx*dx));
			for (int i = Lp*Lp/2; i < Lp + Lp*Lp/2; i ++)
				phi[i] = m + random.nextGaussian()*sqrt((1-m*m)/(dx*dx));
		}
	}
	
	public void readInitialConfiguration(){
		try{
			File myFile = new File("../../../research/javaData/configs/inputConfig");
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
			
			phi[i] = 0.00001*random.nextGaussian()/(dx);
			//phi[i] = m + random.nextGaussian()*sqrt((1-m*m)/(dx*dx));
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
		//double backgroundH = params.fget("H");
		//double stripeH = params.fget("H Stripe");
		//setExternalField(backgroundH, stripeH);
		//dT = params.fget("dT");
		params.set("R/dx", R/dx);
		params.set("Lp", Lp);
		//params.set("Free Energy", freeEnergy);
		
		if(params.sget("Interaction") == "Circle"){
			circleInteraction = true;
		}else{
			circleInteraction = false;
		}
		noiseParameter = params.fget("Noise");
		if(params.sget("Dynamics?") == "Langevin Conserve M")
			magConservation = true;
		else if(params.sget("Dynamics?") == "Langevin No M Conservation")
			magConservation = false;
		theory = params.sget("Approx");
		horizontalSlice = params.fget("Horizontal Slice");
		verticalSlice = params.fget("Vertical Slice");
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
	

	

	
	public void simulateUnstable(){
		convolveWithRange(phi, phi_bar, R);		
		for (int i = 0; i < Lp*Lp; i++){
			delPhi[i] = -dt*(-phi_bar[i]+T* scikit.numerics.Math2.atanh(phi[i])- H);
			phi[i] += delPhi[i];
		}
		t += dt;
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
				dF_dPhi = -phi_bar[i]+T*(phi[i]+pow(phi[i],3)) - H;	
				Lambda[i] = 1;
			}else{
				dF_dPhi = -phi_bar[i]+T* scikit.numerics.Math2.atanh(phi[i])- H;
				//dF_dPhi = -phi_bar[i]+T*(-log(1.0-phi[i])+log(1.0+phi[i]))/2.0 - H;
				if(theory == "HalfStep"){
					Lambda[i] = 1;
					entropy = -((1.0 + phi[i])*log(1.0 + phi[i]) +(1.0 - phi[i])*log(1.0 - phi[i]))/2.0;
				}else{
					Lambda[i] = sqr(1 - phi[i]*phi[i]);				
					entropy = -((1.0 + phi[i])*log(1.0 + phi[i]) +(1.0 - phi[i])*log(1.0 - phi[i]))/2.0;
				}
			}
			delPhi[i] = - dt*Lambda[i]*dF_dPhi + sqrt(Lambda[i]*(dt*2*T)/dx)*noise();
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
	
	public PointSet getHslice(){
		int y = (int) (horizontalSlice * Lp);
		double slice[] = new double[Lp];
		for (int x = 0; x < Lp; x++) {
			slice[x] = phi[Lp*y + x];
		}
		return new PointSet(0, dx, slice);
	}

	public PointSet getVslice(){
		int x = (int) (verticalSlice * Lp);
		double slice[] = new double[Lp];
		for (int y = 0; y < Lp; y++) {
			slice[y] = phi[Lp*y + x];
		}
		return new PointSet(0, dx, slice);
	}

	public PointSet get_delHslice(){
		int y = (int) (horizontalSlice * Lp);
		double slice[] = new double[Lp];
		for (int x = 0; x < Lp; x++) {
			slice[x] = delPhi[Lp*y + x];
		}
		return new PointSet(0, dx, slice);
	}	

	public PointSet get_delVslice(){
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
	
//	private void setExternalField(double backgroundH, double stripeH){
//		for (int i = Lp; i < Lp*Lp; i++)
//			HH[i] = backgroundH;
//		for (int i = 0; i < Lp; i++)
//			HH[i] = stripeH;
//		for (int i = Lp*Lp/2; i < Lp + Lp*Lp/2; i++)
//			HH[i] = stripeH;	
//	}
}

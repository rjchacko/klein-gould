package rachele.ising.dim2;


import scikit.dataset.Accumulator;
import static java.lang.Math.*;
import scikit.numerics.fft.FFT1D;
import scikit.numerics.fft.managed.ComplexDouble2DFFT;

/*
* Calculates the structure factor
*/
public class StructureFactor {
	ComplexDouble2DFFT fft;	// Object to perform transforms
	//RealDoubleFFT_Radix2 fft1D;
	FFT1D fft1d;
	double[] fftData;       // Fourier transform data
	public double sFactor [];
	int Lp;                 // # elements per side
	double L;               // the actual system length, L = Lp*dx, where dx is lattice spacing
	double R;               // characteristic length.  x-axis is k*R.
	double dt;
	double kRmin, kRmax;
	static double squarePeakValue = 4.4934092;
	static double circlePeakValue = 5.135622302;
	int squarePeakInt, circlePeakInt;
	double lastHpeak, lastVpeak, lastCpeak, lastSpeak, lastCmax, lastCmin;
	int noBins = 1024;	
	
	Accumulator accCircle;
	Accumulator accHorizontal;
	Accumulator accVertical;
	Accumulator accAvH;
	Accumulator accAvV;
	Accumulator accAvC;
	Accumulator accPeakH;
	Accumulator accPeakV;
	Accumulator accPeakC;
	Accumulator acc2PeakC;
	Accumulator accPeakHslope;
	Accumulator accPeakVslope;
	Accumulator accPeakCslope;
	Accumulator ringFT;
	Accumulator ringData;
	
	public StructureFactor(int Lp, double L, double R, double kRbinWidth, double dt) {
		this.Lp = Lp;
		this.L = L;
		this.R = R;
		this.dt = dt;
		
		sFactor = new double [Lp*Lp];
		
		kRmin = (2*PI*2/L)*R; // explicitly exclude constant (k=0) mode
		kRmax = (2*PI*(Lp/2)/L)*R;
		
		
		double dblePeakLength = squarePeakValue*L/(2*PI*R);
		squarePeakInt = (int)dblePeakLength;
		//if(abs(2*PI*squarePeakInt*R/L - squarePeakValue) >= abs(2*PI*(squarePeakInt+1)*R/L - squarePeakValue))
			//squarePeakInt = squarePeakInt + 1;
		double kRvalue = R*2*PI*squarePeakInt/L;
		double kRvalue2 = R*2*PI*(squarePeakInt+1)/L;
		System.out.println("square kR = " + kRvalue + " target value = " + squarePeakValue);
		System.out.println("square kR 2 = " + kRvalue2 + " target value = " + squarePeakValue);	
		//squarePeakInt += 1;
		dblePeakLength = circlePeakValue*L/(2*PI*R);
		circlePeakInt = (int)dblePeakLength;
		if(abs(2*PI*circlePeakInt*R/L - circlePeakValue) >= abs(2*PI*(circlePeakInt+1)*R/L - circlePeakValue))
			circlePeakInt = circlePeakInt + 1;
		kRvalue = R*2*PI*circlePeakInt/L;
		System.out.println("circle kR = " + kRvalue + " target value = " + circlePeakValue);		
		System.out.println("circle int = " + circlePeakInt + " Square int = " + squarePeakInt);
			//
//		dblePeakLength = circlePeakValue*L/(2*PI*R);
//		circlePeakInt = (int)dblePeakLength;
//		if(abs(2*PI*circlePeakInt*R/L - circlePeakValue) >= abs(2*PI*(circlePeakInt+1)*R/L - circlePeakValue))
//			circlePeakInt = circlePeakInt + 1;
		
		accCircle = new Accumulator(kRbinWidth);
		accHorizontal = new Accumulator(kRbinWidth);
		accVertical = new Accumulator(kRbinWidth);
		accAvH = new Accumulator(kRbinWidth);
		accAvV = new Accumulator(kRbinWidth);
		accAvC = new Accumulator(kRbinWidth);
		accPeakH = new Accumulator(dt);
		accPeakV = new Accumulator(dt);
		accPeakC = new Accumulator(dt);
		acc2PeakC = new Accumulator(dt);
		accPeakHslope = new Accumulator(dt);
		accPeakVslope = new Accumulator(dt);
		accPeakCslope = new Accumulator(dt);
		ringFT = new Accumulator(.01);
		ringData = new Accumulator(1);
		
		fft = new ComplexDouble2DFFT(Lp, Lp);
		fft1d = new FFT1D(noBins);
		fftData = new double[2*Lp*Lp];
	}

	public double circlekRValue(){
		int circleCount = 0;
		double kRSum = 0;
		for (int y = -Lp/2; y < Lp/2; y++) {
			for (int x = -Lp/2; x < Lp/2; x++) {
				double kR = (2*PI*sqrt(x*x+y*y)/L)*R;
				if(kR >= circlePeakValue - PI*R/(0.5*L) && kR <= circlePeakValue + PI*R/(0.5*L)){
					//System.out.println("Circle kR value = " + kR + "Target value = " + circlePeakValue);
					circleCount += 1;
					kRSum += kR;
				}
			}
		}
		double kRAve = kRSum / circleCount;
		//System.out.println("kR value = " + " " + kRAve);
		return kRAve;
	}
	
	public double getCircleKR(){
		return R*2*PI*circlePeakInt/L;
	}
	
	public Accumulator getPeakV() {
		return accPeakV;
	}
	
	public Accumulator getPeakH() {
		return accPeakH;
	}

	public Accumulator getPeakC() {
		return accPeakC;
	}
	
	public Accumulator get2PeakC() {
		return acc2PeakC;
	}
	
	public Accumulator getPeakVslope() {
		return accPeakCslope;
	}
	
	public Accumulator getPeakHslope() {
		return accPeakCslope;
	}
	
	public Accumulator getPeakCslope() {
		return accPeakCslope;
	}

	public Accumulator getAccumulatorC() {
		return accCircle;
	}

	public Accumulator getAccumulatorH() {
		return accHorizontal;
	}
	
	public Accumulator getAccumulatorV() {
		return accVertical;
	}

	public Accumulator getAccumulatorVA() {
		return accAvV;
	}

	public Accumulator getAccumulatorHA() {
		return accAvH;
	}
	
	public Accumulator getAccumulatorCA() {
		return accAvC;
	}
	
	public double peakValueV(){
		return lastVpeak;
	}
	
	public double peakValueH(){
		return lastHpeak;
	}
	
	public double peakValueC(){
		return lastCpeak;
	}
	
	public double peakValueS(){
		return lastSpeak;
	}
	
	
	public double minC(){
		return lastCmin;
	}
	
	public double maxC(){
		return lastCmax;
	}
	
	
	public double kRmin() {
		return kRmin;
	}
	
	public double kRmax() {
		return kRmax;
	}
	
	public Accumulator getRingFT(){
		return ringFT;
	}

	public Accumulator getRingInput(){
		return ringData;
	}
	
	public void setBounds(double kRmin, double kRmax) {
		this.kRmin = kRmin;
		this.kRmax = kRmax;
	}
	

	
	
	public void accumulateAll(double t, double[] data) {
		double dx = (L/Lp);
		for (int i = 0; i < Lp*Lp; i++) {
			fftData[2*i] = data[i]*dx*dx;
			fftData[2*i+1] = 0;
		}
		accumulateAllAux(t);
	}
	
	public void accumMin(double[] data,double kRcircle){
		
		double dx = (L/Lp);
		for (int i = 0; i < Lp*Lp; i++) {
			fftData[2*i] = data[i]*dx*dx;
			fftData[2*i+1] = 0;
		}
			fft.transform(fftData);
			fftData = fft.toWraparoundOrder(fftData);
			//double dP_dt;
		
	
			for (int i=0; i < Lp*Lp; i++){
				double re = fftData[2*i];
				double im = fftData[2*i+1];
				sFactor[i] = (re*re + im*im)/(L*L);
			}
	
			//Instead of a circular average, we want the structure factor in the vertical and
			//horizontal directions.
			//vertical component
			int y = squarePeakInt;	
			int x=0;
			double kR = (2*PI*sqrt(x*x+y*y)/L)*R;
			int i = Lp*((y+Lp)%Lp) + (x+Lp)%Lp;			
			lastVpeak = sFactor[i];
	
			//horizontal component
			x = squarePeakInt;
			y=0;
			kR = (2*PI*sqrt(x*x+y*y)/L)*R;
			i = Lp*((y+Lp)%Lp) + (x+Lp)%Lp;
			lastHpeak = sFactor[i];
	
			//circularly averaged	
			//double binSize = 2*PI/noBins;
			double [] ringBins = new double [noBins];
			int clumpCt = 0;
			double clumpSum = 0;
			int stripeCt = 0; 
			double stripeSum = 0;
			ringFT.clear();
			for (y = -Lp/2; y < Lp/2; y++) {
				for (x = -Lp/2; x < Lp/2; x++) {
					kR = (2*PI*sqrt(x*x+y*y)/L)*R;
					if (kR >= kRmin && kR <= kRmax) {
						i = Lp*((y+Lp)%Lp) + (x+Lp)%Lp;
						if(kR >= kRcircle - PI*R/(0.50*L) && kR <= kRcircle + PI*R/(0.50*L)){
							double angle = atan((double)abs(y)/(double)abs(x));
							double director = 0.244346095;
								if (x <= 0 && y >=0)
									angle = PI -angle;
								else if (x <=0 && y<=0)
									angle += PI;
								else if (x >= 0 && y <= 0)
									angle = 2*PI - angle;
								if (angle == 2*PI)
									angle = 0;
								double angleS = angle - director;
								if(angleS >= -.1 && angleS <= .1){
									stripeCt +=1;
									stripeSum += sFactor[i];	
									//System.out.println(angleS+ " " + sFactor[i]);
								}else if (angleS >= -.1 + PI && angleS <= .1 + PI){
									stripeCt +=1;
									stripeSum += sFactor[i];	
									//System.out.println(angleS+ " " + sFactor[i]);
								}else if(angleS >= -.1 + PI/3 && angleS <= .1 + PI/3){
									clumpCt += 1;
									clumpSum += sFactor[i];						
									//System.out.println(angleS+ " " + sFactor[i]);
								}else if(angleS >= -.1 + 2*PI/3 && angleS <= .1 + 2*PI/3){
									clumpCt += 1;
									clumpSum += sFactor[i];						
									//System.out.println(angleS+ " " + sFactor[i]);
								}else if(angleS >= -.1 + 4*PI/3 && angleS <= .1 + 4*PI/3){
									clumpCt += 1;
									clumpSum += sFactor[i];						
									//System.out.println(angleS+ " " + sFactor[i]);
								}else if(angleS >= -.1 + 5*PI/3 && angleS <= .1 + 5*PI/3){
									clumpCt += 1;
									clumpSum += sFactor[i];						
									//System.out.println(angleS+ " " + sFactor[i]);
								}
									
								
//								ringFT.accum(angle*360/(2*PI),sFactor[i]);
								ringFT.accum(angle,sFactor[i]);
								//System.out.println(sFactor[i]);
								int j = (int)floor(angle*noBins/(2*PI));
								ringBins[j] += sFactor[i];
//								if(j==52 || j == 564){
//									stripeCt +=1;
//									stripeSum += sFactor[i];
//								}else if(j == 256 || j== 416 || j == 768 || j == 928) {
//									clumpCt += 1;
//									clumpSum += sFactor[i];
//								}

						}
					}	
				}
			}
			lastCpeak = clumpSum/clumpCt;
			lastSpeak = stripeSum/stripeCt;
			ringData.clear();
//			for (int j = 0; j < noBins; j ++){
//				//System.out.println(j + " "+ ringBins[j]);
//				ringData.accum(j,ringBins[j]);
//			}
	}
	
	public double getSF(double [] data){
		double dx = (L/Lp);
		for (int i = 0; i < Lp*Lp; i++) {
			fftData[2*i] = data[i]*dx*dx;
			fftData[2*i+1] = 0;
		}
		fft.transform(fftData);
		fftData = fft.toWraparoundOrder(fftData);
		for (int i=0; i < Lp*Lp; i++){
			double re = fftData[2*i];
			double im = fftData[2*i+1];
			sFactor[i] = (re*re + im*im)/(L*L);
		}
		double sf = 0;
		for (int y = -squarePeakInt; y < 3*squarePeakInt; y = y + 2*squarePeakInt) {
			int x=0;
			int i = Lp*((y+Lp)%Lp) + (x+Lp)%Lp;
			sf += sFactor[i];	
		}		
		for (int x = -squarePeakInt; x < 3*squarePeakInt; x = x + 2*squarePeakInt) {
			int y=0;
			int i = Lp*((y+Lp)%Lp) + (x+Lp)%Lp;
			sf += sFactor[i];	
		}		
		sf /= 4;
		return sf;
	}
	
	public void accumulateAllAux(double t) {
		// compute fourier transform
		fft.transform(fftData);
		fftData = fft.toWraparoundOrder(fftData);
		//double dP_dt;
		
		for (int i=0; i < Lp*Lp; i++){
			double re = fftData[2*i];
			double im = fftData[2*i+1];
			sFactor[i] = (re*re + im*im)/(L*L);
		}
		//Instead of a circular average, we want the structure factor in the vertical and
		//horizontal directions.
		
		//vertical component
		for (int y = -Lp/2; y < Lp/2; y++) {
			int x=0;
			double kR = (2*PI*sqrt(x*x+y*y)/L)*R;
			if (kR >= kRmin && kR <= kRmax) {
				int i = Lp*((y+Lp)%Lp) + (x+Lp)%Lp;
				accVertical.accum(kR, sFactor[i]);
				accAvV.accum(kR, sFactor[i]);
				if(abs(y) == squarePeakInt){
					//accPeakV.accum(t, sFactor[i]);
					//dP_dt = (sFactor[i] - lastVpeak)/dt;
					//lastVpeak = sFactor[i];
					//accPeakVslope.accum(t, dP_dt);
				}
				if(abs(y) == circlePeakInt){
					//accPeakC.accum(t, sFactor[i]);
					//System.out.println("accing");
				}
			}
		}
		
		//horizontal component
		for (int x = -Lp/2; x < Lp/2; x++) {
			int y=0;
			double kR = (2*PI*sqrt(x*x+y*y)/L)*R;
			if (kR >= kRmin && kR <= kRmax) {
				int i = Lp*((y+Lp)%Lp) + (x+Lp)%Lp;
				//double re = fftData[2*i];
				//double im = fftData[2*i+1];
				//double sfValue = (re*re + im*im)/(L*L);
				accHorizontal.accum(kR, sFactor[i]);
				accAvH.accum(kR, sFactor[i]);
				if(abs(x) == squarePeakInt){
					accPeakH.accum(t, sFactor[i]);
					//dP_dt = (sFactor[i] - lastHpeak)/dt;
					lastHpeak = sFactor[i];
					//accPeakHslope.accum(t, dP_dt);
				}
//				if(abs(x) == circlePeakInt){
//					//accPeakC.accum(t, sFactor[i]);
//				}
			}
		}		
	
		//circularly averaged	
		//double [] sfValues = new double [Lp*Lp];
		int count = 0;
		lastCpeak = 0;
		accCircle.clear();
		for (int y = -Lp/2; y < Lp/2; y++) {
			for (int x = -Lp/2; x < Lp/2; x++) {
				double kR = (2*PI*sqrt(x*x+y*y)/L)*R;
				//if (kR >= kRmin && kR <= kRmax) {
				if(kR > 0){
					int i = Lp*((y+Lp)%Lp) + (x+Lp)%Lp;
//					double re = fftData[2*i];
//					double im = fftData[2*i+1];
//					double sfValue = (re*re + im*im)/(L*L);
					//accCircle.accum(kR, sFactor[i]);
					//accAvC.accum(kR, sFactor[i]);
					//sfValues[count] = sFactor[i];

					if(kR >= circlePeakValue - PI*R/(0.5*L) && kR <= circlePeakValue + PI*R/(0.5*L)){
						count += 1;
						//accPeakC.accum(t, sFactor[i]);
						//dP_dt = (sFactor[i] - lastCpeak)/dt;
						lastCpeak += sFactor[i];
						//accPeakCslope.accum(t, dP_dt);
						//System.out.println("Circle kR value = " + kR + "Target value = " + circlePeakValue);
					}
				}
				//}
			}
		}
		lastCpeak /= count;
		sFactor[0] = 0;
		
		//shift sfactor so that origin is in center
		double [] temp = new double [Lp*Lp];
		for (int i = 0; i<Lp*Lp; i++){
			int x = i%Lp;
			int y = i/Lp;
			x += Lp/2; y += Lp/2;
			x = x%Lp; y = y%Lp;
			int j = Lp*((y+Lp)%Lp) + (x+Lp)%Lp;
			temp[j] = sFactor[i];
		}
		for(int i = 0; i<Lp*Lp; i++)
			sFactor[i] = temp[i];
		
	}

	public void takeFT(double[] data){
		double dx = (L/Lp);
		for (int i = 0; i < Lp*Lp; i++) {
			fftData[2*i] = data[i]*dx*dx;
			fftData[2*i+1] = 0;
		}
		fft.transform(fftData);
		fftData = fft.toWraparoundOrder(fftData);
		for (int i=0; i < Lp*Lp; i++){
			double re = fftData[2*i];
			double im = fftData[2*i+1];
			sFactor[i] = (re*re + im*im)/(L*L);
		}		
	}
	
	public void accumulateMelt(boolean circleOn, double[] data, int maxi) {
		takeFT(data);
		if (circleOn==false){
			//vertical component
			//for (int y = -Lp/2; y < Lp/2; y++) {
			//int x=0;
			int x = 0;
			int y = squarePeakInt;
			int i = Lp*((y+Lp)%Lp) + (x+Lp)%Lp;
			//if(abs(y) == squarePeakInt){
			lastVpeak = sFactor[i];
			//}
			//}

			//horizontal component
			x = 1;
			y = squarePeakInt;
			i = Lp*((y+Lp)%Lp) + (x+Lp)%Lp;
			//if(abs(y) == squarePeakInt){
			lastHpeak = sFactor[i];
			
			//			for (int x = -Lp/2; x < Lp/2; x++) {
//				int y=0;
//				int i = Lp*((y+Lp)%Lp) + (x+Lp)%Lp;
//				if(abs(x) == squarePeakInt){
//					lastHpeak = sFactor[i];
//				}
//			}		
		}else{
			//circularly averaged	
			//int count = 0;
			lastCpeak = 0;
			accCircle.clear();
			for (int y = -Lp/2; y < Lp/2; y++) {
				for (int x = -Lp/2; x < Lp/2; x++) {
					int i = Lp*((y+Lp)%Lp) + (x+Lp)%Lp;
					if(i==maxi)
						lastCpeak += sFactor[i];
				}
			}
		}

	}	

	public boolean findStripeDirection(double [] data){
		boolean verticalStripes=true;
		takeFT(data);
			int vert = Lp*((squarePeakInt+Lp)%Lp) + (Lp)%Lp;
			int hor = Lp*((Lp)%Lp) + (squarePeakInt+Lp)%Lp;
			if(sFactor[vert] > sFactor[hor]){
				verticalStripes=false;
			}else{
				verticalStripes=true;
			}
		return verticalStripes;
	}

	public int clumpsOrStripes(double [] data){
		int maxi=1;
		double maxsf=0.0;
		double maxkR = 0.0;
		takeFT(data);
		for (int y = -Lp/2; y < Lp/2; y++) {
			for (int x = -Lp/2; x < Lp/2; x++) {
				double kR = (2*PI*sqrt(x*x+y*y)/L)*R;
				if(kR >= circlePeakValue - PI*R/(0.5*L) && kR <= circlePeakValue + PI*R/(0.5*L)){
				int i = Lp*((y+Lp)%Lp) + (x+Lp)%Lp;
				if(sFactor[i] > maxsf){
					maxsf=sFactor[i];
					maxi=i;
					maxkR=kR;
					System.out.println(sFactor[i]+" "+ maxsf + " maxi " + maxi + " kr " +kR);
				}
				}
			}
		}
		System.out.println(" maxi = " +maxi  + " maxkR = " + maxkR);
		return maxi;
	}
	
	public void accumExact(double t, double[] data, int kR1int, int kR2int) {
		double dx = (L/Lp);
		for (int i = 0; i < Lp*Lp; i++) {
			fftData[2*i] = data[i]*dx*dx;
			fftData[2*i+1] = 0;
		}
		fft.transform(fftData);
		fftData = fft.toWraparoundOrder(fftData);
		//double dP_dt;
		
		for (int i=0; i < Lp*Lp; i++){
			double re = fftData[2*i];
			double im = fftData[2*i+1];
			sFactor[i] = (re*re + im*im)/(L*L);
		}

		//vertical component
		for (int y = -Lp/2; y < Lp/2; y++) {
			int x=0;
			int i = Lp*((y+Lp)%Lp) + (x+Lp)%Lp;
			if(abs(y) == kR1int)	accPeakC.accum(t, sFactor[i]);
			if(abs(y) == kR2int)    acc2PeakC.accum(t, sFactor[i]);
		}

		//horizontal component
		for (int x = -Lp/2; x < Lp/2; x++) {
			int y=0;
			int i = Lp*((y+Lp)%Lp) + (x+Lp)%Lp;
			if(abs(y) == kR1int)	accPeakC.accum(t, sFactor[i]);
			if(abs(y) == kR2int)    acc2PeakC.accum(t, sFactor[i]);
		}		
	}		
}

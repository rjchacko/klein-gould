package kang.ising;

import static java.lang.Math.PI;
import scikit.numerics.fft.managed.ComplexDouble2DFFT;

public class StructureFactorOpt {
	ComplexDouble2DFFT fft;
	double[] fftData;       // Fourier transform data
	public double sFactor [];
	int Lp;                 // # elements per side
	double L;               // the actual system length, L = Lp*dx, where dx is lattice spacing
	static double squarePeakValue = 4.4934092;
	static double circlePeakValue = 5.135622302;

	public StructureFactorOpt(int Lp, double L) {
		this.Lp = Lp;
		this.L = L;
		sFactor = new double [Lp*Lp];	
		fft = new ComplexDouble2DFFT(Lp, Lp);
		fftData = new double[2*Lp*Lp];
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

		shiftSFactor();
	}
	
	public void shiftSFactor(){
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
	
	public double findSquarePeak(double R, boolean xdirection){
		double peakValue;
		int x,y;
		shiftSFactor();
		int kInt = findBestSquareInt(R);
		double kValue = 2*PI*kInt*R/L;
		double difference = kValue - squarePeakValue;
		if(xdirection){
			x = kInt;
			y=0;
			System.out.println("kx = " + kInt + " difference = " + difference);
		}else{
			x=0;
			y = kInt;
			System.out.println("ky = " + kInt + " difference = " + difference);
		}
		int i = Lp*((y+Lp)%Lp) + (x+Lp)%Lp;	
		peakValue = sFactor[i];
		return peakValue;
	}
	
	private int findBestSquareInt(double R){
		double bestValue = L*squarePeakValue/(2.0*PI*R);
		int bestInt = (int)Math.floor(bestValue + 0.5d);
		return bestInt;
	}
}
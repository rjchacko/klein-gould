package rachele.ising.dim2;

import rachele.util.FileUtil;
import scikit.numerics.fft.managed.ComplexDoubleFFT;
import scikit.numerics.fft.managed.ComplexDoubleFFT_Mixed;
import scikit.numerics.fft.managed.ComplexDouble2DFFT;

public class FindCoefficients {
	public double [] a;
	public double [][] aa;
	public double [] inputFunction, phiSlice;
	public int Lp;
	public String fileName;
	ComplexDoubleFFT fft;
	ComplexDouble2DFFT fft2D;
	
	
	public FindCoefficients(int Lp, String fileName){
		this.fileName = fileName;
		this.Lp = Lp;
		a = new double [Lp];
		aa = new double [Lp][Lp];
		fft = new ComplexDoubleFFT_Mixed(Lp);
		fft2D = new ComplexDouble2DFFT(Lp,Lp);
	}
	
	public void findCoefficientsFromFile(){
		phiSlice = new double [Lp];
		phiSlice = FileUtil.readConfigFromFile(fileName, Lp);
		inputFunction = new double [Lp];
		double [] fftScratch = new double [2*Lp];
		//System.out.println("find coeff from file");
		for(int i = 0; i < Lp; i++)
			inputFunction[i] = 1.0/(1-Math.pow(phiSlice[i],2));
		for (int i = 0; i < Lp; i++) {
			fftScratch[2*i] = inputFunction[i];
			fftScratch[2*i+1] = 0;
		}
		fft.transform(fftScratch);
		for (int i = 0; i < Lp; i++) 
			a[i] = fftScratch[2*i];			
	}
	
	public void findCoefficientsFromSlice(double [] slice){
		double [] inputFunction = new double [Lp];
		double [] fftScratch = new double [2*Lp];
		for(int i = 0; i < Lp; i++)
			inputFunction[i] = 1.0/(1-Math.pow(slice[i],2));
		for (int i = 0; i < Lp; i++) {
			fftScratch[2*i] = inputFunction[i];
			fftScratch[2*i+1] = 0;
		}
		fft.transform(fftScratch);
		for (int i = 0; i < Lp; i++) 
			a[i] = fftScratch[2*i];			
	}
	
	public void findCoefficientsFromConfig(double [] config){
		double [] inputFunction = new double [Lp*Lp];
		double [] fftScratch = new double [2*Lp*Lp];
		for(int i = 0; i < Lp*Lp; i++)
			inputFunction[i] = 1.0/(1-Math.pow(config[i],2));
		for (int i = 0; i < Lp*Lp; i++) {
			fftScratch[2*i] = inputFunction[i];
			fftScratch[2*i+1] = 0;
		}
		fft2D.transform(fftScratch);
		for (int i = 0; i < Lp*Lp; i++){
			int kx = i%Lp;
			int ky = i/Lp;
			aa[kx][ky] = fftScratch[2*i];			
		}
	}
	
//	public void findCoefficientsMod(double T, ){
//		double [] phiSlice = new double [Lp];
//		phiSlice = FileUtil.readConfigFromFile(fileName, Lp);
//		double [] G = new double [Lp];
//		double [] F = new double [Lp];
//		for (int i = 0; i < Lp; i++){
//			G[i] = -(1.0-2.0*Math.pow(phiSlice[i],2)+Math.pow(phiSlice[i],4));
//			F[i] = G[i]*T/(1-Math.pow(phiSlice[i],2)) 
//			+ 4 * phiSlice[i]*(1.0-Math.pow(a, phiSlice,2))*();
//		}
//	}

	
}


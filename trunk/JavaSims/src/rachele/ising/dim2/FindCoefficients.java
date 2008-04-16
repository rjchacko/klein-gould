package rachele.ising.dim2;

import rachele.util.FileUtil;
import scikit.dataset.PointSet;
import scikit.numerics.fft.managed.ComplexDoubleFFT;
import scikit.numerics.fft.managed.ComplexDoubleFFT_Mixed;

public class FindCoefficients {
	public double [] a;
	public int Lp;
	private double [] phiSlice, inputFunction, fftScratch;
	public String fileName = "../../../research/javaData/configs1d/config";
	ComplexDoubleFFT fft;
	
	
	public FindCoefficients(int Lp){
		this.Lp = Lp;
		a = new double [Lp];
		inputFunction = new double[Lp];
		phiSlice = new double [Lp];
		fftScratch = new double [2*Lp];
		fft = new ComplexDoubleFFT_Mixed(Lp);
	}
	
	public void findCoefficientsFromFile(){
		phiSlice = FileUtil.readConfigFromFile(fileName, Lp);
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
	
	
	public PointSet getPhiSlice(){
		return new PointSet(0, 1, phiSlice);
	}
}


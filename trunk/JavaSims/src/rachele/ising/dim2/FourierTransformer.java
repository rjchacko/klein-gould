package rachele.ising.dim2;
import scikit.numerics.fft.managed.ComplexDouble2DFFT;
import scikit.numerics.fft.managed.ComplexDoubleFFT;
import scikit.numerics.fft.managed.ComplexDoubleFFT_Mixed;

public class FourierTransformer {

	ComplexDoubleFFT fft1D;
	ComplexDouble2DFFT fft2D;
	public double [] scratch1D;
	public double [] scratch2D;
	public double [] dst1D;
	public double [] dst2D;
	public int L;
	
	
	public FourierTransformer(int L){
		this.L = L;
		fft1D = new ComplexDoubleFFT_Mixed(L);
		fft2D = new ComplexDouble2DFFT(L,L);
		scratch1D = new double[2*L];
		scratch2D = new double[2*L*L];
		dst1D = new double[L];
		dst2D = new double[L*L];
	}

	public double [] calculate1DFT(double [] src){
		for (int i = 0; i < L; i ++){
			scratch1D[2*i] = src[i];
			scratch1D[2*i+1] = 0;
		}
		fft1D.transform(scratch1D);
		for (int i = 0; i < L*L; i ++)
			dst1D[i] = scratch1D[2*i];
		return dst1D;
	}
	
	public double [] calculate2DFT(double [] src){
		for (int i = 0; i < L*L; i ++){
			scratch2D[2*i] = src[i];
			scratch2D[2*i+1] = 0;
		}
		fft2D.transform(scratch2D);
		for (int i = 0; i < L*L; i ++)
			dst2D[i] = scratch2D[2*i];
		return dst2D;
	}
	
	public double [] calculateSF2D(double [] src){
		for (int i = 0; i < L*L; i++) {
			scratch2D[2*i] = src[i];
			scratch2D[2*i+1] = 0;
		}
		fft2D.transform(scratch2D);
		scratch2D = fft2D.toWraparoundOrder(scratch2D);
	
		for (int i=0; i < L*L; i++){
			double re = scratch2D[2*i];
			double im = scratch2D[2*i+1];
			dst2D[i] = (re*re + im*im)/(L*L);
		}
		return dst2D;
	}
	
}

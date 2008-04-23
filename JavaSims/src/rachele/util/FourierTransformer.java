package rachele.util;
import scikit.numerics.fft.managed.ComplexDouble2DFFT;
import scikit.numerics.fft.managed.ComplexDoubleFFT;
import scikit.numerics.fft.managed.ComplexDoubleFFT_Mixed;
import scikit.numerics.fn.Function2D;

public class FourierTransformer {

	ComplexDoubleFFT fft1D;
	ComplexDouble2DFFT fft2D;
	public double [] scratch1D;
	public double [] scratch2D;
	public double [] scratch2D2;
	public double [] dst1D;
	public double [] dst2D;
	public int L;
	
	
	public FourierTransformer(int L){
		this.L = L;
		fft1D = new ComplexDoubleFFT_Mixed(L);
		fft2D = new ComplexDouble2DFFT(L,L);
		scratch1D = new double[2*L];
		scratch2D = new double[2*L*L];
		scratch2D2 = new double[2*L*L];
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
	
	public double [] calculateSF2D(double [] src, boolean centered, boolean zeroCenter){
		dst2D = find2DSF(src);
		if (zeroCenter) dst2D[0] = 0;
		if (centered) center(dst2D);
		return dst2D;
	}
	
	public void center(double [] src){
		double [] temp = new double [L*L];
		for (int i = 0; i<L*L; i++){
			int x = i%L;
			int y = i/L;
			x += L/2; y += L/2;
			x = x%L; y = y%L;
			int j = L*((y+L)%L) + (x+L)%L;
			temp[j] = src[i];
		}
		for(int i = 0; i<L*L; i++)
			src[i] = temp[i];	
	}
	
	public double [] convolve1DwithFunction(double [] src, Function2D fn){
		for (int i = 0; i < L; i++){
			scratch1D[2*i] = src[i];
			scratch1D[2*i + 1] = 0;
		}
		fft1D.transform(scratch1D);
		scratch1D = fft1D.toWraparoundOrder(scratch1D);		
		for (int i = 0; i < L; i++){
			scratch1D[2*i] *= -1*fn.eval((double)i,1.0);
			scratch1D[2*i+1] *= -1*fn.eval((double)i,1.0);
		}
		fft1D.backtransform(scratch1D);
		for (int i = 0; i < L; i++){
			dst1D[i] = scratch1D[2*i]/(L*L);
			//System.out.println(dst1D[i]);
		}
		return dst1D;
	}
	
	
	public double [] convolve2D(double [] src1, double [] src2){
		
		for (int i = 0; i < L*L; i++) {
			scratch2D[2*i] = src1[i];
			scratch2D[2*i+1] = 0;
			scratch2D2[2*i] = src2[i];
			scratch2D2[2*i+1] = 0;
		}		
		fft2D.transform(scratch2D);
		fft2D.transform(scratch2D2);
		scratch2D = fft2D.toWraparoundOrder(scratch2D);
		scratch2D2 = fft2D.toWraparoundOrder(scratch2D2);
		
		for (int i = 0; i < 2 * L * L; i++)
			scratch2D[i] *= scratch2D2[i];
		
		fft2D.backtransform(scratch2D);
		for (int i = 0; i < L * L; i++)		
			dst2D[i] = scratch2D[2*i]/(L*L);
		return dst2D;
	}
	
	private double [] find2DSF(double [] src){
		double [] dst = new double [L*L];
		for (int i = 0; i < L*L; i++) {
			scratch2D[2*i] = src[i];
			scratch2D[2*i+1] = 0;
		}
		fft2D.transform(scratch2D);
		scratch2D = fft2D.toWraparoundOrder(scratch2D);
	
		for (int i=0; i < L*L; i++){
			double re = scratch2D[2*i];
			double im = scratch2D[2*i+1];
			dst[i] = (re*re + im*im)/(L*L);
		}
		return dst;
	}
}	



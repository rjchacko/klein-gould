package rachele.util;
import scikit.numerics.fft.managed.ComplexDouble2DFFT;
import scikit.numerics.fft.managed.ComplexDoubleFFT;
import scikit.numerics.fft.managed.ComplexDoubleFFT_Mixed;
import scikit.numerics.fn.Function1D;
import scikit.numerics.fn.Function2D;

public class FourierTransformer {

	ComplexDoubleFFT fft1D;
	ComplexDouble2DFFT fft2D;
	public double [] scratch1D;
	public double [] scratch1D2;
	public double [] scratch2D;
	public double [] scratch2D2;
	public int L;
	
	
	public FourierTransformer(int L){
		this.L = L;
		fft1D = new ComplexDoubleFFT_Mixed(L);
		fft2D = new ComplexDouble2DFFT(L,L);
		scratch1D = new double[2*L];
		scratch1D2 = new double[2*L];
		scratch2D = new double[2*L*L];
		scratch2D2 = new double[2*L*L];
	}

	public double [] calculate1DFT(double [] src){
		double [] dst = new double[L];
		for (int i = 0; i < L; i ++){
			scratch1D[2*i] = src[i];
			scratch1D[2*i+1] = 0;
		}
		fft1D.transform(scratch1D);
		for (int i = 0; i < L; i ++)
			dst[i] = scratch1D[2*i]/L;
		return dst;
	}

	public double [] calculate1DBackFT(double [] src){
		double [] dst = new double[L];
		for (int i = 0; i < L; i ++){
			scratch1D[2*i] = src[i];
			scratch1D[2*i+1] = 0;
		}
		fft1D.transform(scratch1D);
		for (int i = 0; i < L; i ++)
			dst[i] = scratch1D[2*i]/L;
		return dst;
	}
	
	public double [] calculate2DFT(double [] src){
		double [] dst = new double [L*L];
		for (int i = 0; i < L*L; i ++){
			scratch2D[2*i] = src[i];
			scratch2D[2*i+1] = 0;
		}
		fft2D.transform(scratch2D);
		for (int i = 0; i < L*L; i ++)
			dst[i] = scratch2D[2*i];
		return dst;
	}
	
	public double [] calculate2DSF(double [] src, boolean centered, boolean zeroCenter){
		double [] dst = new double [L*L];
		dst = find2DSF(src);
		if (zeroCenter) dst[0] = 0;
		if (centered) center(dst);
		return dst;
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
	
	public double [] convolve1DwithFunction(double [] src, Function1D fn){
		double [] dst = new double [L];
		for (int i = 0; i < L; i++){
			scratch1D[2*i] = src[i];
			scratch1D[2*i + 1] = 0;
		}
		fft1D.transform(scratch1D);
		for (int x = -L/2; x < L/2; x++) {
			int i = (x+L)%L;
			scratch1D[2*i] *=  fn.eval((double)x);
			scratch1D[2*i+1] *=  fn.eval((double)x);
		}
		fft1D.backtransform(scratch1D);
		for (int i = 0; i < L; i++)
			dst[i] = scratch1D[2*i]/(L);
		return dst;
	}

	public double [] backConvolve1DwithFunction(double [] src, Function1D fn){
		double [] dst = new double [L];
		for (int i = 0; i < L; i++){
			scratch1D[2*i] = src[i];
			scratch1D[2*i + 1] = 0;
		}
		fft1D.backtransform(scratch1D);
		for (int x = -L/2; x < L/2; x++) {
			int i = (x+L)%L;
			scratch1D[2*i] *=  fn.eval((double)x);
			scratch1D[2*i+1] *=  fn.eval((double)x);
		}
		fft1D.transform(scratch1D);
		for (int i = 0; i < L; i++)
			dst[i] = scratch1D[2*i]/(L);
		return dst;
	}
	
	public double [] backConvolve1D(double [] src1, double [] src2){
		double [] dst = new double [L];
		src1 = calculate1DBackFT(src1);
		src2 = calculate1DBackFT(src2);
		for (int i = 0; i < L; i++)
			dst[i] = src1[i]*src2[i];
		dst = calculate1DFT(dst);
		return dst;		
	}
	
	public double [] convolve1D(double [] src1, double [] src2){
		double [] dst = new double [L];
		src1 = calculate1DFT(src1);
		src2 = calculate1DFT(src2);
		for (int i = 0; i < L; i++)
			dst[i] = src1[i]*src2[i];
		dst = calculate1DBackFT(dst);
		return dst;	
	}
	
	public double [] convolve2D(double [] src1, double [] src2){
		double [] dst = new double [L*L];
		for (int i = 0; i < L*L; i++) {
			scratch2D[2*i] = src1[i];
			scratch2D[2*i+1] = 0;
			scratch2D2[2*i] = src2[i];
			scratch2D2[2*i+1] = 0;
		}		
		fft2D.backtransform(scratch2D);
		fft2D.backtransform(scratch2D2);
		scratch2D = fft2D.toWraparoundOrder(scratch2D);
		scratch2D2 = fft2D.toWraparoundOrder(scratch2D2);
		
		for (int i = 0; i < 2 * L * L; i++)
			scratch2D[i] *= scratch2D2[i];
		
		fft2D.backtransform(scratch2D);
		for (int i = 0; i < L * L; i++)		
			dst[i] = scratch2D[2*i]/(L*L);
		return dst;
	}

	public double [] convolve2DwithFunction(double [] src, Function2D fn){
		double [] dst = new double [L*L];
		for (int i = 0; i < L*L; i++){
			scratch2D[2*i] = src[i];
			scratch2D[2*i + 1] = 0;
		}
		fft2D.transform(scratch2D);
		for (int y = -L/2; y < L/2; y++) {
			for (int x = -L/2; x < L/2; x++) {
				int i = (x+L)%L + L*((y+L)%L);
				scratch2D[2*i] *=  fn.eval((double)x, (double)y);
				scratch2D[2*i+1] *=  fn.eval((double)x, (double)y);
			}
		}
		fft2D.backtransform(scratch2D);
		for (int i = 0; i < L*L; i++)
			dst[i] = scratch2D[2*i]/(L*L);
		return dst;
	}
	
	public double [] find2DSF(double [] src, double size){
		double dx = size/L;
		double [] dst = new double [L*L];
		for (int i = 0; i < L*L; i++) {
			scratch2D[2*i] = src[i]*dx*dx;
			scratch2D[2*i+1] = 0;
		}
		fft2D.transform(scratch2D);
		scratch2D = fft2D.toWraparoundOrder(scratch2D);
	
		for (int i=0; i < L*L; i++){
			double re = scratch2D[2*i];
			double im = scratch2D[2*i+1];
			dst[i] = (re*re + im*im)/(size*size);
		}
		return dst;		
		
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



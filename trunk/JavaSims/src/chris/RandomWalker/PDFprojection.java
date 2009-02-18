package chris.RandomWalker;

import scikit.jobs.Control;
import scikit.jobs.Simulation;
import scikit.jobs.params.DirectoryValue;
import scikit.jobs.params.StringValue;

public class PDFprojection extends Simulation {

	private double lambda, dx, x[], pdf[];
	private int dim;
	private Pdf2step pfn;
	
	public static void main(String[] args) {
		new Control(new PDFprojection(), "Parameters");
	}
	
	public PDFprojection(){
		
		PDFprojectionConstructor();
		return;
	}
	
	public void PDFprojectionConstructor(){
		
		
		return;
	}
	
	public void animate() {
		
	}

	public void clear() {
		
	}

	public void load(Control c) {
		params.add("Data Directory",new DirectoryValue("/Users/cserino/Desktop/"));
		params.add("File Name", new StringValue("foo"));
		params.add("Dimension", (int) 2);
		params.add("dx", 1.01e-4);
		params.add("lambda", 0.5);
	}

	public void run() {
	
		double gsum;
		int Npts;
		double lnow;
		
		dim    = params.iget("Dimension");
		lambda = params.fget("lambda");
		dx     = params.fget("dx");
		gsum   = 1/(1-lambda);
		Npts   = (int)(2*gsum/dx);
		Npts   = Npts + 1 - (Npts % 2);
		x      = new double[Npts];
		pdf    = new double[Npts];
		lnow   = lambda;
		
		pfn.setDim(dim);
		
		x[(int)(Npts/2)+1]   = 0;
		pdf[(int)(Npts/2)+1] = pfn.eval(x[(int)(Npts/2)+1]);
		
		for (int jj = 0 ; jj < (int)(Npts/2) ; jj++){
			x[jj]                   = dx*jj - gsum;
			x[jj+(int)(Npts/2)+1]   = dx*(jj+1);
			pdf[jj]                 = pfn.eval(x[jj]);
			pdf[jj+(int)(Npts/2)+1] = pfn.eval(x[jj+(int)(Npts/2)+1]);
		}
		
		
		while(lnow > dx){
			
			// convolve pdf w/ pfn(x*lnow) 
			// and store result in pdf
			
			lnow = lambda*lnow;
		}
	
		
		// print x and pdf(x) to file
	
		return;
	}
	

/*	public void convolve(double[] src, double[] dst, Function1D fn) {
		double L1 = dim1*dx1;
		for (int i = dim1-1; i >= 0; i--) {
			scratch[2*i+0] = src[i]*dx1;
			scratch[2*i+1] = 0;
		}
		
		fft.transform(scratch);
		scratch = fft.toWraparoundOrder(scratch);

		for (int x = -dim1/2; x < dim1/2; x++) {
			int i = (x+dim1)%dim1;
			double k = 2*PI*x/L1;
			double J = fn.eval(k);
			scratch[2*i+0] *= J;
			scratch[2*i+1] *= J;
		}
		
		fft.backtransform(scratch);
		for (int i = 0; i < dim1; i++) {
			dst[i] = scratch[2*i+0] / L1;
		}
	}*/
	
	
	
}

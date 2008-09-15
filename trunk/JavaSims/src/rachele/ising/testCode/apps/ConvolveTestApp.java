package rachele.ising.testCode.apps;

import java.awt.Color;

import scikit.dataset.PointSet;
import scikit.graphics.dim2.Grid;
import scikit.graphics.dim2.Plot;
import scikit.jobs.Control;
import scikit.jobs.Job;
import scikit.jobs.Simulation;
import rachele.util.FourierTransformer;

public class ConvolveTestApp extends Simulation{


	Plot initial = new Plot("init");
	Plot finalPlot = new Plot("final");
	Plot fourierSpace = new Plot("fourier");
	Plot convolution = new Plot("convolution");
	Plot direct = new Plot("Direct");
	Plot convolve = new Plot("Convolve");
	Grid grid = new Grid("grid");
	Grid grid2 = new Grid("grid2");
	Grid grid3 = new Grid("grid3");
	int size = 128;
	FourierTransformer ft = new FourierTransformer(size);
	double [] iA = new double [size];
	double [] iB = new double [size];
	double [] fA = new double [size];
	double [] fB = new double [size];
	double [] iAB = new double [size];
	double [] kA = new double [size];
	double [] kB = new double [size];
	double [] kAB = new double [size];
	double [] conkAkB = new double [size];
	double [] finalConvolve = new double [size];
	double [] finalDirect = new double [size];
	double [] config2d = new double [size*size];
	double [] config2dB = new double [size*size];
	double [] config2di = new double [size*size];
	
	double [] config2dk= new double [size*size];
	double [] config2df = new double [size*size];
	double [] config2dconvolvef= new double [size*size];
	
	public static void main(String[] args) {
		new Control(new ConvolveTestApp(), "Convolution Test");
	}
	
	public void animate() {
		
		initial.registerLines("init A", new PointSet(0, 1, iA), Color.BLACK);
//		initial.registerLines("init B", new PointSet(0, 1, iB), Color.BLACK);
		fourierSpace.registerLines("FT A", new PointSet(0, 1, kA), Color.BLUE);
//		fourierSpace.registerLines("FT B", new PointSet(0, 1, kB), Color.BLUE);
		finalPlot.registerLines("final A", new PointSet(0, 1, fA), Color.BLACK);
//		finalPlot.registerLines("final B", new PointSet(0, 1, fB), Color.BLUE);
//		initial.registerLines("init AB", new PointSet(0,1,iAB), Color.RED);
	
//		fourierSpace.registerLines("FT of AB", new PointSet(0,1,kAB), Color.BLACK);
		convolution.registerLines("convolution of AB", new PointSet(0,1,conkAkB), Color.BLUE);
		
//		finalPlot.registerLines("final convolution", new PointSet(0,1,finalConvolve), Color.BLUE);
//		finalPlot.registerLines("final direct", new PointSet(0,1,finalDirect), Color.BLACK);
		
		direct.registerLines("init AB", new PointSet(0,1,iAB), Color.BLACK);
		direct.registerLines("FT of AB", new PointSet(0,1,kAB), Color.BLUE);
		direct.registerLines("final direct", new PointSet(0,1,finalDirect), Color.RED);
		
		convolve.registerLines("init AB", new PointSet(0,1,iAB), Color.BLACK);
		convolve.registerLines("convolution of AB", new PointSet(0,1,conkAkB), Color.BLUE);
		convolve.registerLines("final convolution", new PointSet(0,1,finalConvolve), Color.RED);
		
		grid.registerData(size, size, config2df);
		grid2.registerData(size, size, config2di);
		grid3.registerData(size, size, config2dconvolvef);
	
	}
	
	


	public void clear() {		
	}

	public void load(Control c) {
		c.frameTogether("plots", initial, finalPlot, fourierSpace, convolution);
		c.frameTogether("compare", convolve, direct, grid, grid2, grid3);
	}

	public void run() {

		double TA = (double)size/2.0;
		double TB = (double)size/3.0;
		double TC = (double)size/4.0;
		double TD = (double)size/5.0;
		for (int i = 0; i < size; i++){
			iA[i] = Math.sin(2*Math.PI*i/TA)*Math.sin(2*Math.PI*i/TB);
			iB[i] = Math.sin(2*Math.PI*i/TC)*Math.sin(2*Math.PI*i/TD);
		}
		for (int i = 0; i < size*size; i++){
			config2d [i] = iA[i%size]*iB[i/size];
			config2dB [i] = iA[i%size]*iA[i/size];
		}
		int index = size*size-1;
		System.out.println("init " + config2d[index] + " fin " + config2df[index]);		
		
		
		// 1D FT then backFT works to give back original function
		kA = ft.calculate1DFT(iA);
		kB = ft.calculate1DFT(iB);
		fA = ft.calculate1DBackFT(kA);
		fB = ft.calculate1DBackFT(kB);
		
		for (int i = 0; i < size; i++)
			iAB[i] = iA[i]*iB[i];
		kAB = ft.calculate1DFT(iAB);
		finalDirect = ft.calculate1DBackFT(kAB);
		
//		regular convolve doesn't work.  must back convolve
//		conkAkB = ft.convolve1D(kA, kB); 
		conkAkB = ft.backConvolve1D(kA, kB);
		finalConvolve = ft.calculate1DBackFT(conkAkB);

		// 2D FT then backFT works to give back original function
		config2dk = ft.calculate2DFT(config2d);
		config2df = ft.calculate2DBackFT(config2dk);
		System.out.println("init " + config2d[index] + " fin " + config2df[index]);	
		
		for (int i = 0; i < size*size; i++){
			config2di [i] = config2d[i]*config2dB[i];
		}
		//convolve doesn't work- must backconvolve
		config2df = ft.calculate2DFT(ft.calculate2DBackFT(config2di));
		config2dconvolvef = ft.calculate2DBackFT(ft.backConvolve2D(ft.calculate2DFT(config2d), ft.calculate2DFT(config2dB)));
		System.out.println("init " + config2di[index] + " fin " + config2dconvolvef[index]);	
		Job.animate();
		
		while(true){
			
		}
		
	}

}

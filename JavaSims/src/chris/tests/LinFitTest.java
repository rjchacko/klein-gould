package chris.tests;

import java.util.Random;

import scikit.dataset.Histogram;
import scikit.graphics.dim2.Plot;
import scikit.jobs.Control;
import scikit.jobs.Job;
import scikit.jobs.Simulation;
import chris.util.FitUtil;
import chris.util.PrintUtil;

public class LinFitTest extends Simulation{

	
	@SuppressWarnings("unused")
	private double fooX[], fooY[], m, b, dm, db, chis, fitY[], width;
	private int N;
	private Random rand;
	private FitUtil fitter;
	@SuppressWarnings("unused")
	private Histogram hData, hFit;
	private Plot graphD = new Plot("Data");
	private Plot graphF = new Plot("Fit");
	
	
	public static void main(String[] args) {
		new Control(new LinFitTest(), "Test LinFit");
	}
	
	public void animate() {
//		graphF.registerLines("Data",hFit,Color.BLUE);
//		graphD.registerPoints("Data",hData,Color.RED);
//		record();
		return;
	}

	public void clear() {
		graphD.clear();
		graphF.clear();
		return;
	}

	public void load(Control c) {
		//params.add("Input Directory",new DirectoryValue("/Users/cserino/Desktop/"));
		c.frameTogether("Fit Test",graphD, graphF);
		return;
	}

	public void run() {
			
		
		PrintUtil.printlnToFile("/Users/cserino/Desktop/TestLinFit.txt","slope","d(slope)","intercept", "d(intercept)", "chi^2", "width");
		rand   = new Random(0);
		fitter = new FitUtil();
		
		double m0, b0, chi2;
		
		while(true){
		
			N = 10000;
			
			fooX  = new double[N];
			fooY  = new double[N];
			fitY  = new double[N];
			hData = new Histogram(0.5);
			hFit  = new Histogram(0.5);
			width = 10;
			
			m0 = 2;
			b0 = 20;
			m  = 0;
			b  = 0;
			
			
//			for (int jj = 0 ; jj < N ; jj++){
//				fooX[jj] = jj;
//				fooY[jj] = m0*jj + b0 + width*rand.nextGaussian();
//			}

			generateDataFile(m0,b0,rand,"/Users/cserino/Desktop/FitData.txt");
			

			fitter.fit("/Users/cserino/Desktop/FitData.txt",1,1,2,width);
			
			m    = fitter.getSlope();
			b    = fitter.getIntercept();
			chi2 = fitter.getChis();
			
//			for (int jj = 0 ; jj < N ; jj++){
//				fitY[jj] = m*fooX[jj]+b;
//				hData.accum(fooX[jj],fooY[jj]);
//				hFit.accum(fooX[jj],fitY[jj]);
//			}
			
			System.out.print("m ="+"\t");
			System.out.println(m);
			System.out.print("b ="+"\t");
			System.out.println(b);
			System.out.print("chis / N ="+"\t");
			System.out.println(chi2/N);		
			
			Job.animate();
			
		}
	}
	
	private void generateDataFile(double m0, double b0, Random rand, String fout){
		
		PrintUtil.printlnToFile(fout,"x","y");
		for (int jj = 0 ; jj < N ; jj++){
			PrintUtil.printlnToFile(fout,jj,m0*jj + b0 + width*rand.nextGaussian());
		}
	}
	
	public void record(){
		
		PrintUtil.printlnToFile("/Users/cserino/Desktop/TestLinFit.txt",m,dm,b, db, chis, width);

	}

}

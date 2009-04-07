package chris.tests;

import java.awt.Color;
import java.util.Random;

import chris.util.PrintUtil;
import chris.util.FitUtil;

import scikit.dataset.Histogram;
import scikit.graphics.dim2.Plot;
import scikit.jobs.Control;
import scikit.jobs.Job;
import scikit.jobs.Simulation;

public class LinFitTest extends Simulation{

	
	private double fooX[], fooY[], m, b, dm, db, chis, fitY[], width;
	private int N;
	private Random rand;
	private FitUtil fitter;
	private Histogram hData, hFit;
	private Plot graphD = new Plot("Data");
	private Plot graphF = new Plot("Fit");
	
	public static void main(String[] args) {
		new Control(new LinFitTest(), "Test LinFit");
	}
	
	public void animate() {
		graphF.registerLines("Data",hFit,Color.BLUE);
		graphD.registerPoints("Data",hData,Color.RED);
		record();
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
		rand = new Random(0);
		
		while(true){
		
			N = 200;
			
			fooX  = new double[N];
			fooY  = new double[N];
			fitY  = new double[N];
			hData = new Histogram(0.5);
			hFit  = new Histogram(0.5);
			
			fitter = new FitUtil(N);
			width = 50;

			for (int jj = 0 ; jj < N ; jj++){
				fooX[jj] = jj;
				fooY[jj] = width*rand.nextGaussian() + jj + 12;
			}

			double temp[] = fitter.fit(fooX, fooY, width,false);
			
			m    = temp[0];
			b    = temp[1];
			dm   = temp[2];
			db   = temp[3];
			chis = temp[4];
			
			for (int jj = 0 ; jj < N ; jj++){
				fitY[jj] = m*fooX[jj]+b;
				hData.accum(fooX[jj],fooY[jj]);
				hFit.accum(fooX[jj],fitY[jj]);
			}
			
			Job.animate();
			
		}
	}
	
	public void record(){
		
		PrintUtil.printlnToFile("/Users/cserino/Desktop/TestLinFit.txt",m,dm,b, db, chis, width);

	}

}

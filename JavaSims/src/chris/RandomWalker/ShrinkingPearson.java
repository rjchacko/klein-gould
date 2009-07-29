package chris.RandomWalker;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;

import scikit.jobs.params.Parameters;
import chris.util.Random;


public class ShrinkingPearson {

	private int pdfR[], pdfX[], bins;
	private double dr, gsum;
	private Random rand;
	
	public ShrinkingPearson(Parameters params){
		
		ShrinkingPearson_Constructor(params);
		return;
	}
	
	public void ShrinkingPearson_Constructor(Parameters params){
		
		double lambda;
		int seed;
		
		lambda = params.fget("lambda");
		dr     = params.fget("Resolution");
		seed   = params.iget("Random Seed");
		
		rand = new Random(seed);
		gsum = 1 / (1 - lambda);
		bins = (int)(Math.ceil(gsum/dr));
		pdfR = new int[bins];
		pdfX = new int[2*bins];
	}
	
	
	public void nextWalk(double lambda, double steps, int d){
		
		double[] r = new double[d];
		double[] rs;
		double tmp = 0;
		
		for (int jj = 0 ; jj < steps ; jj++){
			rs = rand.nextSpherePt(d,Math.pow(lambda,jj));
			for (int kk = 0 ; kk < d ; kk++){
				r[kk] += rs[kk];
				tmp += rs[kk]*rs[kk];
			}
			tmp = 0;
		}
		
		add2pdf(r, d);
		return;
	}
	
	
	private void add2pdf(double[] r, int d){
		double rr = 0;
		
		for (int jj = 0 ; jj < d ; jj++){
			rr += r[jj]*r[jj];
		}
		rr = Math.sqrt(rr);
		
		pdfR[(int)(rr/dr)]++;
		pdfX[(int)((r[0]+gsum)/dr)]++;
		return;
	}

	public void printPDF(String fout, int steps, int nw, double lambda){
	
		try{
			File file = new File(fout);
			PrintWriter pw = new PrintWriter(new FileWriter(file, true), true);
			pw.print("%% lambda = ");
			pw.print(lambda);
			pw.print("   %% n_walkers = ");
			pw.print(nw);
			pw.print("   %% n_steps = ");
			pw.println(steps);
			pw.print("r");
			pw.print("\t");
			pw.print("P(r)");
			pw.print("\t");
			pw.print("dP(r)");
			pw.print("\t");
			pw.print("x_j");
			pw.print("\t");
			pw.print("P(x_j)");
			pw.print("\t");
			pw.print("dP(x_j)");
			pw.println();
			for (int jj = 0 ; jj < bins ; jj++){
				pw.print(dr*jj);
				pw.print("\t");
				pw.print((double)(pdfR[jj])/nw);
				pw.print("\t");
				pw.print(Math.sqrt(pdfR[jj])/nw);
				pw.print("\t");
				pw.print(dr*jj-gsum);
				pw.print("\t");
				pw.print((double)(pdfX[jj])/nw);
				pw.print("\t");
				pw.print(Math.sqrt(pdfX[jj])/nw);
				pw.println();
			}			
			for (int jj = bins ; jj < 2*bins ; jj++){
				pw.print(-1);
				pw.print("\t");
				pw.print(-1);
				pw.print("\t");
				pw.print(-1);
				pw.print("\t");
				pw.print(dr*jj-gsum);
				pw.print("\t");
				pw.print((double)(pdfX[jj])/nw);
				pw.print("\t");
				pw.print(Math.sqrt(pdfX[jj])/nw);
				pw.println();
			}	
			pw.close();
		}
		catch (IOException ex){
			ex.printStackTrace();
		}
		return;
	}
	
}

package chris.RandomWalker;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import chris.util.CopyUtil;
import chris.util.MathUtil;
import chris.util.PrintUtil;
import scikit.jobs.Control;
import scikit.jobs.Simulation;
import scikit.jobs.params.ChoiceValue;
import scikit.jobs.params.DirectoryValue;
import scikit.jobs.params.DoubleValue;
import scikit.jobs.params.Parameters;
import scikit.jobs.params.StringValue;

public class ProbProp2D extends Simulation{
	
	private double dx, os, pdfold[], pdfnew[], gsum, ss[], lambda;
	private int Npts, Niter, jmin, jmax;
	private String fout;

	public static void main(String[] args){
		new Control(new ProbProp2D(), "Parameters");		
	}

	public void load(Control c) {
		
		params.add("Data Directory",new DirectoryValue("/Users/cserino/Desktop"));
		params.add("File Name", new StringValue("default"));
		params.add("Calculate In", new ChoiceValue("Real Space","Fourier Space"));
		params.add("Resolution", (double) 0.01);
		params.add("lambda", new DoubleValue(0.5, 0., 1.));
		params.add("Status");
	}

	public void run() {
		
		setup(params);
		
		for (int jj = 0 ; jj < Niter ; jj++){
			params.set("Status",(int)(100*jj/Niter));
			getss(jj);
			for (int kk = 0 ; kk < Npts ; kk++){
				if(pdfold[kk] == 0) continue;
				for (int ll = kk-jmin; ll <= kk+jmax; ll++){
					pdfnew[ll] = pdfold[kk]*ss[ll];
				}
			}
			pdfold = CopyUtil.copyArray(pdfnew, Npts);
			pdfnew = new double[Npts];
		}
		
		printData();
		// print to file
	}
	
	private void setup(Parameters params){
		
		gsum   = params.fget("lambda");
		gsum   = 1/(1-gsum);
		dx     = params.fget("Resolution");
		Npts   = (int)(Math.ceil(2*gsum/dx));
		Npts   = Npts + (1-Npts%2); 	// make the number of points odd 
									//so zero is at the center of a bin
		os     = (Npts + dx)/2;
		lambda = params.fget("lambda");
		Niter  = (int)(Math.log(dx)/Math.log(lambda));
		pdfold = new double[Npts];
		pdfnew = new double[Npts];
		fout   = params.sget("Data Directory") + File.separator + params.sget("File Name") + ".txt";
	
		pdfold[(int)(Npts/2)] = 1;	// walk begins at origin with unit probability
		PrintUtil.printlnToFile(fout,"x","pdf");
		return;		
	}
	
	private void getss(int p){
		
		double lp = MathUtil.pow(lambda, p);
		double norm = 0;
		int jmin, jmax;
		jmin  = (int)((-lp+os)/dx);
		jmax  = (int)(( lp+os)/dx);
		if((jmax-jmin)%2==0){
			ss = null; // by symmetry, width should be odd
			return;
		}
		
		for (int jj = jmin ; jj <= jmax ; jj++){
			ss[jj-jmin] = Math.pow((1-MathUtil.pow(getx(jj),2)),-0.5);
			norm        += ss[jj-jmin];
		}
		for (int jj = jmin ; jj <= jmax ; jj++){
			ss[jj-jmin] = ss[jj-jmin]/norm;
		}
		
		return;
	}
	
	private double getx(int index){
		
		return index*dx - os;
	}
	
	private void printData(){

		try{
			File file = new File(fout);
			PrintWriter pw = new PrintWriter(new FileWriter(file, true), true);
			for (int jj = 0 ; jj < Npts ; jj++){
				pw.print(jj*dx+os);
				pw.print("\t");
				pw.print(pdfold[jj]);
				pw.println();
			}			
			pw.close();
		}
		catch (IOException ex){
			ex.printStackTrace();
		}
		return;
	}

	public void animate() {
		
	}

	public void clear() {
		
	}
	
}

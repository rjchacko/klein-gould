package chris.RandomWalker;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import chris.util.CopyUtil;
import chris.util.MathUtil;
import chris.util.PrintUtil;
import scikit.jobs.Control;
import scikit.jobs.Job;
import scikit.jobs.Simulation;
import scikit.jobs.params.ChoiceValue;
import scikit.jobs.params.DirectoryValue;
import scikit.jobs.params.DoubleValue;
import scikit.jobs.params.Parameters;
import scikit.jobs.params.StringValue;

public class ProbProp2D extends Simulation{
	
	private double dx, dx2, pdfold[], pdfnew[], gsum, step[], lambda;
	private int Npts, mpt, Nsteps, jmax;
	private String fout, pout;

	public static void main(String[] args){
		new Control(new ProbProp2D(), "Parameters");		
	}

	public void load(Control c) {
		
		params.add("Data Directory",new DirectoryValue("/Users/cserino/Documents/BU/Research/Pearson"));
		params.add("File Name", new StringValue("default"));
		params.add("Calculate In", new ChoiceValue("Real Space","Fourier Space"));
		params.add("Resolution", (double) 1e-4);
		params.add("lambda", new DoubleValue(0.5, 0., 1.));
		params.add("Status");
	}

	public void run() {
		
		setup(params);
		firststep();

		
		for (int ii = 2 ; ii <= Nsteps ; ii++){
			getstep(ii);
			// debug
			System.out.println(ii);
			// end debug
			params.set("Status",(int)(100*ii/Nsteps));
			Job.animate();
			for (int jj = 0 ; jj < Npts ; jj++){
				if(pdfold[jj] == 0) continue;
//				// debug
//				double pa = 0;
//				// end debug
				for (int kk = (-jmax + 1) ; kk < jmax ; kk++){
					pdfnew[jj+kk] += pdfold[jj]*step[sgn(kk)];
//					// debug
//					pa += pdfold[jj]*step[kk];
//					// end debug
				}
//				// debug
//				PrintUtil.printlnToFile(pout,pa-pdfold[jj]);
//				// end debug

			}
			pdfold = CopyUtil.copyArray(pdfnew, Npts);
			pdfnew = new double[Npts];
			// debug
			double norm = 0;
			for(int mm = 0 ; mm < Npts ; mm++){
				norm += pdfold[mm];
			}
			PrintUtil.printlnToFile(pout,-ii,norm);
			// end debug
		}
		
		printData();
		params.set("Status","Done");
		Job.signalStop();
		Job.animate();	
		
		return;
	}
	
	private void setup(Parameters params){

		params.set("Status","Setup");

		lambda = params.fget("lambda");
		gsum   = 1/(1-lambda);
		dx     = params.fget("Resolution");
		dx2    = dx/2-dx/2000; 	// takes the point at the center of the bin minus
								// a small offset so that the step -/-> infinity 
		Npts   = (int)(2*Math.ceil(gsum/dx));
		Npts   = Npts + (Npts+1)%2; // makes Npts odd ==> symmetry about zero
		mpt    = (int)(Npts/2);
		Nsteps = (int)(Math.log(dx)/Math.log(lambda));
		pdfold = new double[Npts];
		pdfnew = new double[Npts];
		fout   = params.sget("Data Directory") + File.separator + params.sget("File Name") + ".txt";
		pout   = params.sget("Data Directory") + File.separator + "Params_" + params.sget("File Name") + ".txt";
			
		PrintUtil.printlnToFile(fout,"x","pdf");
		PrintUtil.printlnToFile(pout,params.toString());
		PrintUtil.printlnToFile(pout,"Nstep = ",Nsteps);
		return;		
	}
	
	private void firststep(){

		params.set("Status","0");
		
		getstep(1);

		for(int jj = (-jmax+1) ; jj < jmax ; jj++){
			pdfold[mpt+jj] = step[sgn(jj)];
		}
		
		return;
	}
	
	private void getstep(int sn){
		
		double lp, norm;
		
		norm = 0;
		if(sn > 1){
			lp = MathUtil.pow(lambda, sn-1);
		}
		else{
			lp = 1;
		}
		
		jmax = (int)(lp/dx - 0.5) + 1;
		if(sn == 1) step = new double[jmax];
		
		for(int jj = 0 ; jj < jmax ; jj++){
			step[jj] = Math.pow(1-MathUtil.pow(getx(jj)/lp, 2),-0.5);
			norm     += step[jj];
			// debugging
			if(step[jj] >= 0){
				if(step[jj] == Double.POSITIVE_INFINITY) PrintUtil.printlnToFile(pout,step[jj]);

			}
			else{
				PrintUtil.printlnToFile(pout,step[jj]);
				PrintUtil.printlnToFile(pout,sn);
				PrintUtil.printlnToFile(pout,lp);
				PrintUtil.printlnToFile(pout,jj);
				PrintUtil.printlnToFile(pout,getx(jj));
				PrintUtil.printlnToFile(pout,getx(jj)/lp);
				PrintUtil.printlnToFile(pout,MathUtil.pow(getx(jj)/lp,2));
				PrintUtil.printlnToFile(pout,1-MathUtil.pow(getx(jj)/lp,2));
				PrintUtil.printlnToFile(pout,Math.pow(1-MathUtil.pow(getx(jj)/lp, 2),-0.5));
			}
			// end debugging
		}
		// this is only the norm for x > 0 so . . . 
		norm = norm - step[0];
		norm = 2*norm;
		norm = norm + step[0];
		
		for(int jj = 0 ; jj < jmax ; jj++){
			step[jj] = step[jj]/(norm);
		}
		// debugging
		norm = 0;
		for(int jj = (-jmax+1) ; jj < jmax ; jj++){
			norm += step[sgn(jj)];
		}
		PrintUtil.printlnToFile(pout,sn,norm);
		// end debugging

		return;
	}
	
	private double getx(int index){
		
		return (index == mpt) ? 0 : index*dx + dx2;
	}
	
	private int sgn(int val){
		
		return (val >= 0) ? val : -val; 
	}
	
	private void printData(){

		double normchk = 0;
		
		try{
			File file = new File(fout);
			PrintWriter pw = new PrintWriter(new FileWriter(file, true), true);
			for (int jj = 0 ; jj < Npts ; jj++){
				pw.print((jj-mpt)*dx);
				pw.print("\t");
				pw.print(pdfold[jj]);
				normchk += pdfold[jj];
				pw.println();
			}			
			pw.close();
		}
		catch (IOException ex){
			ex.printStackTrace();
		}
		PrintUtil.printlnToFile(pout,"Final Normalization = ", normchk);
		return;
	}

	public void animate() {
		
		return;
	}

	public void clear() {
		
		return;
	}
	
}

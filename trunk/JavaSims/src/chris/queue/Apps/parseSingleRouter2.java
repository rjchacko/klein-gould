package chris.queue.Apps;

import java.io.BufferedInputStream;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.FilenameFilter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.text.DecimalFormat;

import scikit.dataset.DatasetBuffer;
import scikit.dataset.Histogram;
import scikit.jobs.Control;
import scikit.jobs.Job;
import scikit.jobs.Simulation;
import scikit.jobs.params.DirectoryValue;

public class parseSingleRouter2 extends Simulation{

	DecimalFormat fmtsd = new DecimalFormat("00.0");
	DecimalFormat fmts  = new DecimalFormat("000");

	int tl, L;
	
	public static void main(String[] args) {
		new Control(new parseSingleRouter2(), "Merge Single Router Data");
	}
	
	public void load(Control c) {
		
		params.add("singleRouterApp Directory", new DirectoryValue("/Users/cserino/Desktop/"));
		params.add("l",(int)5);
		params.add("\u03BB", 0.19);
		params.add("Time Series Length", (int)500000);
		params.add("Status");

	}
	
	public void animate() {
		
		return;
	}

	public void clear() {

		return;
	}

	public void run() {
		
		int NN;
		Histogram hN, hnj[];
		
		tl            = params.iget("Time Series Length");
		L             = params.iget("l");
	    hN            = new Histogram(1);
		hnj           = new Histogram[L]; 	
		double[][] n  = new double[L][tl];
		double[][] n2 = new double[L][tl];
		String pth    = params.sget("singleRouterApp Directory");
		String tpth[] = new String[]{pth+"/ts",pth+"/nt",pth+"/nj",pth+"/njn"};
		
		for(int jj = 0 ; jj < L ; jj++){
			hnj[jj]  = new Histogram(1);
		}

		
		// Generate a list of files 
		
		// 1. time series data
		File dirContents = new File(tpth[0]);
		String[] dta     = dirContents.list();

		// 2. N histogram data
		dirContents  = new File(tpth[1]);
		String shN[] = dirContents.list();
		
		// 3. nj histogram data
		String shnj[][] = new String[L][];
		dirContents    = new File(tpth[2]);
		for(int jj = 0 ; jj < L ; jj++){
			final String tmp      = fmts.format(jj);
			FilenameFilter filter = new FilenameFilter() {
				public boolean accept(File dir, String name) {
					return name.contains("_histn"+tmp); 
				}
			};
			shnj[jj] = dirContents.list(filter);
		}
		
		// Read in data 

		NN = Math.min(dta.length, shN.length);
		
		FileInputStream fis;
		BufferedInputStream bis;
		BufferedReader bir;
		String rin;
		int pd;
		for (int jj = 0 ; jj < NN ; jj++){
		
			// 1. time series data
			try {
				fis  = new FileInputStream(tpth[0]+File.separator+dta[jj]);
				bis  = new BufferedInputStream(fis);
				bir  = new BufferedReader(new InputStreamReader(bis));
				rin  = bir.readLine(); // skip header
				int count = 0;
				while ( (rin = bir.readLine()) != null ){
					n[count%L][count/L]  += Integer.parseInt(rin);	
					n2[count%L][count/L] += n[count%L][count/L]*n[count%L][count/L];	
					count++;
				}
			} catch (FileNotFoundException e) {
				e.printStackTrace();
			} catch (IOException e) {
				e.printStackTrace();
			}		
		
			// 2. N histogram data
			try {
				fis = new FileInputStream(tpth[1]+File.separator+shN[jj]);
				bis = new BufferedInputStream(fis);
				bir = new BufferedReader(new InputStreamReader(bis));
				while ( (rin = bir.readLine()) != null ){
					pd = rin.indexOf('\t');
					hN.accum(Double.parseDouble(rin.substring(0,pd)),Double.parseDouble(rin.substring(pd+1)));
				}
			} catch (FileNotFoundException e) {
				e.printStackTrace();
			} catch (IOException e) {
				e.printStackTrace();
			}
		
			// 3. nj histogram data
			for (int kk = 0 ; kk < L ; kk++){
				try {
					fis = new FileInputStream(tpth[2]+File.separator+shnj[kk][jj]);
					bis = new BufferedInputStream(fis);
					bir = new BufferedReader(new InputStreamReader(bis));
					while ( (rin = bir.readLine()) != null ){
						pd = rin.indexOf('\t');
						hnj[kk].accum(Double.parseDouble(rin.substring(0,pd)),Double.parseDouble(rin.substring(pd+1)));
					}
				} catch (FileNotFoundException e) {
					e.printStackTrace();
				} catch (IOException e) {
					e.printStackTrace();
				}
			}

			if(jj % 5 == 0){
				params.set("Status",fmtsd.format(100.*jj/NN)+"% done");
				Job.animate();
			}
		}
		
		params.set("Status","Writing data");
		Job.animate();

		// 1. time series data
		File file;
		PrintWriter pw;
		try{
			file = new File(pth+File.separator+"TimeSeries.txt");
			pw = new PrintWriter(new FileWriter(file, true), true);
			pw.println("t"+"\t"+"<n>_0"+"\t"+"var(n_0)"+"\t"+"<n>_1"+"\t"+"var(n_1)"+"\t"+" .  .  .  ");
			for (int jj = 0 ; jj < tl ; jj++){
				pw.print(jj);
				pw.print("\t");
				for (int kk = 0 ; kk < L ; kk++){
					n[kk][jj]  /= NN;
					n2[kk][jj] = (n2[kk][jj] - NN*n[kk][jj]*n[kk][jj])/NN;
					// now n[j][t]  is the average  over the NN trials of n_j(t) where 0 <= j < l
					// now n2[j][t] is the variance over the NN trials of n_j(t) where 0 <= j < l
					pw.print(n[kk][jj]);
					pw.print("\t");
					pw.print(n2[kk][jj]);
					pw.print("\t");
				}
				pw.println();
			}
			for (int jj = 0 ; jj < tl ; jj++){
				pw.println(n[jj]);
			}			
			pw.close();
		}
		catch (IOException ex){
			ex.printStackTrace();
		}
		
		// 2. N histogram data
		DatasetBuffer hh;
		try{
			hh = hN.copyData();
			file = new File(pth+File.separator+"Nhist.txt");
			pw = new PrintWriter(new FileWriter(file, true), true);
			pw.print("N_trials = ");
			pw.println(NN);
			for (int jj = 0 ; jj < hh.size(); jj++){
				pw.print(hh.x(jj));
				pw.print("\t");
				pw.println(hh.y(jj));
			}			
			pw.close();
		}
		catch (IOException ex){
			ex.printStackTrace();
		}
		
		// 3. nj histogram data
		for (int ii = 0 ; ii < L ; ii++){
			try{
				hh = hnj[ii].copyData();
				file = new File(pth+File.separator+"n"+fmts.format(ii)+"hist.txt");
				pw = new PrintWriter(new FileWriter(file, true), true);
				pw.print("N_trials = ");
				pw.println(NN);
				for (int jj = 0 ; jj < hh.size(); jj++){
					pw.print(hh.x(jj));
					pw.print("\t");
					pw.println(hh.y(jj));
				}			
				pw.close();
			}
			catch (IOException ex){
				ex.printStackTrace();
			}
		}
		
		params.set("Status","Done");
		Job.signalStop();
		Job.animate();
		return;
	}
}
	

	
	
	
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

public class parseSingleRouter3 extends Simulation{

	DecimalFormat fmtsd = new DecimalFormat("00.0");
	DecimalFormat fmts  = new DecimalFormat("000");
	DecimalFormat fmtL  = new DecimalFormat("0000000");	

	int L;
	
	public static void main(String[] args) {
		new Control(new parseSingleRouter3(), "Merge Single Router Data");
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
		
		L               = params.iget("l");
		String pth      = params.sget("singleRouterApp Directory");
		String tpth     = pth+"/njn";
		int Nmax        = (int)(3*L/(1-params.fget("\u03BB")*L));
		int Norm[]      = new int[Nmax];
		Histogram h[][] = new Histogram[L][Nmax];
		String s[][][]  = new String[L][Nmax][];
		
		for (int jj = 0 ; jj < Nmax ; jj++)
			for (int kk = 0 ; kk < L ; kk++)
				h[kk][jj] = new Histogram(1);
		
		
		// make list of files
		FilenameFilter filter;
		File dirContent = new File(tpth);
		for(int jj = 0 ; jj < Nmax ; jj++){
			for(int kk = 0 ; kk < L ; kk++){
				final String tmp = fmts.format(kk)+"g"+fmts.format(jj);
				filter = new FilenameFilter() {
					public boolean accept(File dir, String name) {
						return name.contains(tmp); 
					}
				};
			s[kk][jj] = dirContent.list(filter);
			}
		}
		
		// get norms
		for(int jj = 0 ; jj < Nmax ; jj++){
				final String tmp = "g"+fmts.format(jj);
				filter = new FilenameFilter() {
					public boolean accept(File dir, String name) {
						return name.contains(tmp); 
					}
				};
			Norm[jj] = dirContent.list(filter).length;
		}
		
		// read in data
		FileInputStream fis;
		BufferedInputStream bis;
		BufferedReader bir;
		String rin;
		int pd;
		for(int jj = 0 ; jj < Nmax ; jj++){
			for(int kk = 0 ; kk < L ; kk++){
				for(int ll = 0 ; ll < s[kk][jj].length ; ll++){
					try {
						fis = new FileInputStream(tpth+File.separator+s[kk][jj][ll]);
						bis = new BufferedInputStream(fis);
						bir = new BufferedReader(new InputStreamReader(bis));
						while ( (rin = bir.readLine()) != null ){
							pd = rin.indexOf('\t');
							h[kk][jj].accum(Double.parseDouble(rin.substring(0,pd)),Double.parseDouble(rin.substring(pd+1)));
						}
					} catch (FileNotFoundException e) {
						e.printStackTrace();
					} catch (IOException e) {
						e.printStackTrace();
					}
				}
			}
			if(jj % 5 == 0){
				params.set("Status",fmtsd.format(100.*jj/Nmax)+"% done");
				Job.animate();
			}
		}
	
		// write data.
		params.set("Status","Writing data");
		Job.animate();
		DatasetBuffer hh;
		File file;
		PrintWriter pw;
		for(int jj = 0 ; jj < Nmax ; jj++){
			for(int kk = 0 ; kk < L ; kk++){
				try{
					hh = h[kk][jj].copyData();
					file = new File(pth+File.separator+"n"+fmts.format(kk)+"g"+fmts.format(jj)+"hist.txt");
					pw = new PrintWriter(new FileWriter(file, true), true);
					pw.print("N_trials = ");
					pw.println(Norm[jj]);
					for (int ll = 0 ; ll < hh.size(); ll++){
						pw.print(hh.x(ll));
						pw.print("\t");
						pw.println(hh.y(ll));
					}			
					pw.close();
				}
				catch (IOException ex){
					ex.printStackTrace();
				}
			}
		}	
		params.set("Status","Done");
		Job.signalStop();
		Job.animate();
		return;
	}	
}
	
	
	
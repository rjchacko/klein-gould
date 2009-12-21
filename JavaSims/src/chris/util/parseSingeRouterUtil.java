package chris.util;

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

public class parseSingeRouterUtil extends Simulation{

	DecimalFormat fmtS = new DecimalFormat("00.0");
	int tl;
	
	public static void main(String[] args) {
		new Control(new parseSingeRouterUtil(), "Merge Single Router Data");
	}
	
	public void load(Control c) {
		
		params.add("singleRouterApp Directory", new DirectoryValue("/Users/cserino/Desktop/"));
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
		
		tl          = params.iget("Time Series Length");
		Histogram h = new Histogram(1);
		int[] n     = new int[tl];
		String pth  = params.sget("singleRouterApp Directory");
		int NN;
		
		// Generate a list of files 
		
		File dirContents = new File(pth);
		    
		FilenameFilter txtfilter = new FilenameFilter() {
			public boolean accept(File dir, String name) {
				return name.endsWith(".txt") && !name.contains("_hist");
			}
		};

		FilenameFilter hstfilter = new FilenameFilter() {
			public boolean accept(File dir, String name) {
				return name.endsWith(".txt") && name.contains("_hist");
			}
		};
		
		String[] dta = dirContents.list(txtfilter);
		String[] hst = dirContents.list(hstfilter);

		NN = dta.length;
		
		for (int jj = 0 ; jj < NN ; jj++){
			try {
				FileInputStream fis = new FileInputStream(pth+File.separator+dta[jj]);
				BufferedInputStream bis = new BufferedInputStream(fis);
				BufferedReader bir = new BufferedReader(new InputStreamReader(bis));
				String rin;
				rin = bir.readLine(); // skip header
				int count = 0;
				while ( (rin = bir.readLine()) != null ){
					n[count++] += Integer.parseInt(rin);
				}
			} catch (FileNotFoundException e) {
				e.printStackTrace();
			} catch (IOException e) {
				e.printStackTrace();
			}			
			try {
				FileInputStream fis = new FileInputStream(pth+File.separator+hst[jj]);
				BufferedInputStream bis = new BufferedInputStream(fis);
				BufferedReader bir = new BufferedReader(new InputStreamReader(bis));
				String rin;
				int pd;
				while ( (rin = bir.readLine()) != null ){
					pd = rin.indexOf('\t');
					h.accum(Double.parseDouble(rin.substring(0,pd)),Double.parseDouble(rin.substring(pd+1)));
				}
			} catch (FileNotFoundException e) {
				e.printStackTrace();
			} catch (IOException e) {
				e.printStackTrace();
			}
			if(jj % 10 == 0){
				params.set("Status",fmtS.format(100.*jj/NN)+"% done");
				Job.animate();
			}
		}
			
		params.set("Status","Writing data");
		Job.animate();

		try{
			File file = new File(pth+File.separator+"TimeSeries.txt");
			PrintWriter pw = new PrintWriter(new FileWriter(file, true), true);
			for (int jj = 0 ; jj < tl ; jj++){
				pw.println(n[jj]);
			}			
			pw.close();
		}
		catch (IOException ex){
			ex.printStackTrace();
		}
		try{
			DatasetBuffer hh = h.copyData();
			File file = new File(pth+File.separator+"Histogram.txt");
			PrintWriter pw = new PrintWriter(new FileWriter(file, true), true);
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
		
		params.set("Status","Done");
		Job.signalStop();
		Job.animate();
		return;
	}
}
	

	
	
	
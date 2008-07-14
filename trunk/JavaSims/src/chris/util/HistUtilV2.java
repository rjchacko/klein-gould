package chris.util;

import java.awt.Color;
import java.io.BufferedInputStream;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStreamReader;

import scikit.dataset.Histogram;
import scikit.graphics.dim2.Plot;
import scikit.jobs.Control;
import scikit.jobs.Job;
import scikit.jobs.Simulation;
import scikit.jobs.params.ChoiceValue;
import scikit.jobs.params.FileValue;

public class HistUtilV2 extends Simulation{

	Histogram histGLOB;
	Plot grid = new Plot("Histogram");
	boolean loglog;
	
	public static void main(String[] args) {
		new Control(new HistUtilV2(), "Create 1D Histogram from TAB");
	}
	
	public void load(Control c) {
		
		params.add("TAB File", new FileValue("/Users/cserino/Desktop/"));
		params.add("Number of Header Rows / Rows to Skip", (int) 1);
		params.add("Column Number", (int) 1);
		params.add("Bin Width", (double) 1);
		params.add("Log-Log", new ChoiceValue("Yes","No"));
		params.add("Status");
		
		c.frame(grid);
		
	}
	
	public void animate() {
		if(loglog) grid.setLogScale(true,true);
		grid.registerPoints("Histogram",histGLOB,Color.RED);
	}

	public void clear() {
		grid.clear();
	}

	public void run() {
		
		params.set("Status","Initializing . . .");
		
		if(params.sget("Log-Log").equals("Yes")){
			loglog = true;
		}
		else{
			loglog = false;
		}

		int cn = params.iget("Column Number");
		int skip = params.iget("Number of Header Rows / Rows to Skip");
		double bw = params.fget("Bin Width");
		
		String finSTRING = params.sget("TAB File");
		File fin = new File(finSTRING);		
		String fout = setupOutFile(finSTRING);

		params.set("Status","Fetching Data . . .");
		double[] data = ReadIn(fin, cn, skip);
		
		params.set("Status","Filling Histogram . . .");
		double[] hist = fillHist(data,bw);
	
		if(loglog){
			params.set("Status","Writing Histogram . . .");
			PrintUtil.printlnToFile(fout,"Coordinate","Number","Log[coord]","Log[#]");
			for (int jj = 0 ; jj < hist.length ; jj=jj+2){
				double coord = hist[jj];
				double num = hist[jj+1];
				PrintUtil.printlnToFile(fout,coord,num,Math.log(coord),Math.log(num));
			}
		}
		else{
			params.set("Status","Writing Histogram . . .");
			PrintUtil.printlnToFile(fout,"Coordinate","Number");
			for (int jj = 0 ; jj < hist.length ; jj=jj+2){
				PrintUtil.printlnToFile(fout,hist[jj],hist[jj+1]);
			}
		}
		Job.animate();
		
		params.set("Status","Done");
	}
	
	private String setupOutFile(String fin){

		int findslash = fin.indexOf('/');
		String ret = fin.substring(0,findslash+1);
		fin = fin.substring(findslash+1);
		findslash = fin.indexOf('/');
		while(findslash != -1){
			ret+=fin.substring(0,findslash+1);
			fin = fin.substring(findslash+1);
			findslash = fin.indexOf('/');
		}
		findslash = fin.indexOf('.');
		ret+=fin.substring(0,findslash);
		ret+="_hist.txt";
		
		return ret;
	}
	
	private double[] fillHist(double[] data, double bw){
		histGLOB = new Histogram(bw);
		for (int jj = 0 ; jj < data.length ; jj++){
			histGLOB.accum(data[jj]);
		}
		
		return histGLOB.copyData();
		
		
	}
	
	private double[] ReadIn(File fin, int cn, int skip){
		
		double[] values = new double[2000000];
		int counter = 0;
		
		try {
			
			FileInputStream fis = new FileInputStream(fin);
			BufferedInputStream bis = new BufferedInputStream(fis);
			BufferedReader bir = new BufferedReader(new InputStreamReader(bis));
			
			String rin;
			int pd;
			
			for (int jj = 0 ; jj < skip ; jj++){
				rin = bir.readLine();
			}
			
			while ( (rin = bir.readLine()) != null ){
				pd = rin.indexOf('\t');
				for(int jj = 1 ; jj < cn ; jj++){
					rin = rin.substring(pd + 1);
					pd = rin.indexOf('\t');
				}
				if(pd == -1){
					values[counter++] = Double.parseDouble(rin);
				}
				else{
					values[counter++] = Double.parseDouble(rin.substring(0,pd));	
				}
			}
			
			params.set("Status","Done");
			
		} catch (FileNotFoundException e) {
			e.printStackTrace();
			params.set("Status","Error!");
		} catch (IOException e) {
			e.printStackTrace();
			params.set("Status","Error!");
		}
		
		double[] ret = new double[counter];
		
		for (int jj = 0 ; jj < counter ; jj++){
			ret[jj] = values[jj];
		}
		
		return ret;

	}
	

}

package chris.util;

import java.awt.Color;
import java.io.BufferedInputStream;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStreamReader;

import scikit.dataset.DatasetBuffer;
import scikit.dataset.Histogram;
import scikit.graphics.dim2.Plot;
import scikit.jobs.Control;
import scikit.jobs.Job;
import scikit.jobs.Simulation;
import scikit.jobs.params.ChoiceValue;
import scikit.jobs.params.FileValue;

public class HistUtil extends Simulation{

	private Histogram histGLOB;
	private DatasetBuffer hist;
	private Plot grid = new Plot("Histogram");
	private boolean loglog, semilog, done;
	private static int dlength = 2000000;
	
	public static void main(String[] args) {
		new Control(new HistUtil(), "Create 1D Histogram from TAB");
	}
	
	public void load(Control c) {
		
		params.add("TAB File", new FileValue("/Users/cserino/Desktop/"));
		params.add("Number of Header Rows / Rows to Skip", (int) 1);
		params.add("Column Number", (int) 1);
		params.add("Bin Width", (double) 1);
		params.add("Scaling Form", new ChoiceValue("Power Law","Exponential"));
		params.add("Status");
		
		c.frame(grid);
		
	}
	
	public void animate() {
		
		if(!done) return;
		if(loglog) grid.setLogScale(true,true);
		if(semilog) grid.setLogScale(false,true);
		grid.registerPoints("Histogram",histGLOB,Color.RED);
		return;
	}

	public void clear() {
		
		grid.clear();
		return;
	}

	public void run() {
		
		while(true){
		
			done = false;
		
			params.set("Status","Initializing . . .");

			if(params.sget("Scaling Form").equals("Power Law")){
				loglog  = true;
				semilog = false;
			}
			else if(params.sget("Scaling Form").equals("Exponential")){
				loglog  = false;
				semilog = true;
			}
			else{

			}

			int cn    = params.iget("Column Number");
			int skip  = params.iget("Number of Header Rows / Rows to Skip");
			double bw = params.fget("Bin Width");
			histGLOB  = new Histogram(bw);

			String finSTRING = params.sget("TAB File");
			File fin = new File(finSTRING);		
			String fout = setupOutFile(finSTRING);

			params.set("Status","Fetching and Binning Data . . .");
			Job.animate();
			ReadIn(fin, cn, skip, bw);
			
			params.set("Status","Writing Histogram . . .");
			Job.animate();
			hist = histGLOB.copyData();
			if(loglog){
				PrintUtil.printlnToFile(fout,"Coordinate","Number","Log[coord]","Log[#]");
				for (int jj = 0 ; jj < hist.size() ; jj++){
					double coord = hist.x(jj);
					double num = hist.y(jj);
					PrintUtil.printlnToFile(fout,coord,num,Math.log(coord),Math.log(num));
				}
			}
			else if (semilog){
				PrintUtil.printlnToFile(fout,"Coordinate","Number","Log[#]");
				for (int jj = 0 ; jj < hist.size() ; jj++){
					double coord = hist.x(jj);
					double num = hist.y(jj);
					PrintUtil.printlnToFile(fout,coord,num,Math.log(num));
				}
			}
			else{
				PrintUtil.printlnToFile(fout,"Coordinate","Number");
				for (int jj = 0 ; jj < hist.size() ; jj++){
					PrintUtil.printlnToFile(fout,hist.x(jj),hist.y(jj));
				}
			}
			
			done = true;
			params.set("Status","Done");
			Job.animate();	
			
			Job.manualStop();
			
		}
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
	
	private void fillHist(double[] data, double bw){
		
		for (int jj = 0 ; jj < data.length ; jj++){
			histGLOB.accum(data[jj]);
		}
		return;	
	}
	
	private void ReadIn(File fin, int cn, int skip, double bw){
		
		double[] values = new double[dlength];
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
				if(counter%dlength == 0){
					fillHist(values,bw);
					counter = 0;
				}
			}
			if(counter != 0){
				double[] foo = CopyUtil.copyArray(values,counter);
				fillHist(foo,bw);
			}
		} 
		catch (FileNotFoundException e) {
			e.printStackTrace();
			params.set("Status","Error!");
		} 
		catch (IOException e) {
			e.printStackTrace();
			params.set("Status","Error!");
		}
		
		return;

	}
	

}

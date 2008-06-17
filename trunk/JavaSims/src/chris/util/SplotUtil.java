package chris.util;

import java.io.BufferedInputStream;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintWriter;

import scikit.dataset.Histogram;
import scikit.jobs.Control;
import scikit.jobs.Simulation;
import scikit.jobs.params.FileValue;

public class SplotUtil extends Simulation{

	Histogram histGLOB;
	
	public static void main(String[] args) {
		new Control(new SplotUtil(), "Create for GNUPLOT Splot from TAB");
	}
	
	public void load(Control c) {
		
		params.add("TAB File", new FileValue("/Users/cserino/Desktop/"));
		params.add("Number of Header Rows / Rows to Skip", (int) 1);
		params.add("Status");
				
	}
	
	public void animate() {
	}

	public void clear() {
	}

	public void run() {
		
		params.set("Status","Initializing . . .");

		int skip = params.iget("Number of Header Rows / Rows to Skip");
		
		String finSTRING = params.sget("TAB File");
		File fin = new File(finSTRING);		
		String fout = setupOutFile(finSTRING);

		params.set("Status","Running . . .");
		
		try {
			FileInputStream fis;
			fis = new FileInputStream(fin);
			BufferedInputStream bis = new BufferedInputStream(fis);
			BufferedReader bir = new BufferedReader(new InputStreamReader(bis));

			String rin;
			int pd;

			int mastercounter = 0;

			while((rin = bir.readLine()) != null){
				
				if(mastercounter++ < skip) continue;
				
				double[] values = new double[1000000];
				int counter = 0;


				while((pd = rin.indexOf('\t')) != -1){
					values[counter++] = Double.parseDouble(rin.substring(0,pd));
					rin = rin.substring(pd+1);
				}
				
				//System.out.println(rin);
				//values[counter++] = Double.parseDouble(rin);
				double[] pout = CopyArray.copyArray(values, counter);
				
				printRow(fout,pout);

			}
		}
		catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

		params.set("Status","Done");
		return;
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
		ret+="_splot.txt";
		
		return ret;
	}
	
	private void printRow(String fout, double[] vals){
		
		double time = vals[0];
			
		try{
			File file = new File(fout);
			PrintWriter pw = new PrintWriter(new FileWriter(file, true), true);
			
			
			for (int jj = 1 ; jj < vals.length; jj++){
				pw.print(time);
				pw.print("\t");
				pw.print(jj);
				pw.print("\t");
				pw.print(vals[jj]);
				pw.println();
			}
			
			pw.println();
		}
		catch (IOException ex){
			ex.printStackTrace();
		}
		
		return;

	}
	

}

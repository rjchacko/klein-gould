package chris.util;

import java.text.DecimalFormat;

import scikit.jobs.Control;
import scikit.jobs.Job;
import scikit.jobs.Simulation;
import scikit.jobs.params.FileValue;

public class ManySplotUtil extends Simulation{

	DecimalFormat fmt = new DecimalFormat("0.00");
	
	public static void main(String[] args) {
		new Control(new ManySplotUtil(), "Create for GNUPLOT Splot from TAB");
	}
	
	public void load(Control c) {
		
		params.add("TAB File", new FileValue("/Users/cserino/Desktop/"));
		params.add("Number of Header Rows / Rows to Skip", (int) 0);
		params.add("Column Number for Y Data", (int) 2);
		params.add("Column Number for Z Data", (int) 3);
		params.add("Status");
				
	}
	
	public void animate() {
	}

	public void clear() {
	}

	public void run() {
		
		
		String bsnm = getBasename();
		String fout = bsnm + "_splot.txt";
		
		for(double ii = 0.5 ; ii < 1 ; ii+= 0.01){
			
			params.set("Status",fmt.format(ii));
			Job.animate();
			
			String fin = bsnm+ "_" + fmt.format(ii) + ".txt";
			
			ReadInUtil rin = new ReadInUtil(fin);
			double[][] tmp = rin.getData(new int[] {1,2,3},params.iget("Number of Header Rows / Rows to Skip"));			
			PrintUtil.printArray4Splot(fout, tmp, 1, 2, tmp[0].length, ii);
			
		}
		
//		ReadInUtil rin = new ReadInUtil(params.sget("TAB File"));
//		double[][] tmp = rin.getData(new int[] {1,2,3},params.iget("Number of Header Rows / Rows to Skip"));
//		String tmpfout = "/Users/cserino/Desktop/TESTOUT.txt";
//		PrintUtil.printArrayToFile(tmpfout,tmp,tmp[0].length,tmp.length);
		
		return;
	}

	private String getBasename(){
		
		String fin = params.sget("TAB File");
		
		int findslash = fin.indexOf('/');
		String ret = fin.substring(0,findslash+1);
		fin = fin.substring(findslash+1);
		findslash = fin.indexOf('/');
		while(findslash != -1){
			ret+=fin.substring(0,findslash+1);
			fin = fin.substring(findslash+1);
			findslash = fin.indexOf('/');
		}
		findslash = fin.indexOf('_');
		ret+=fin.substring(0,findslash);
		
		return ret;

	}

}

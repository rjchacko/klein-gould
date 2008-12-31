package chris.RandomWalker;

import java.text.DecimalFormat;

import chris.util.PrintUtil;
import chris.util.ReadInUtil;

import scikit.jobs.Control;
import scikit.jobs.Job;
import scikit.jobs.Simulation;
import scikit.jobs.params.FileValue;

public class GeometricRsqr extends Simulation{

	DecimalFormat fmt = new DecimalFormat("0.00");
	
	public static void main(String[] args) {
		new Control(new GeometricRsqr(), "get <r^2> from distribution");
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
		String fout = bsnm + "_R2.txt";
		PrintUtil.printlnToFile(fout,"lambda","<R^2>");
		
		for(double ii = 0.5 ; ii < 1 ; ii+= 0.01){
			
			params.set("Status",fmt.format(ii));
			Job.animate();
			
			String fin = bsnm+ "_" + fmt.format(ii) + ".txt";
			
			ReadInUtil rin = new ReadInUtil(fin);
			double[][] tmp = rin.getData(new int[] {1,2,3},params.iget("Number of Header Rows / Rows to Skip"));			
			double R2 = getR2(tmp);
			PrintUtil.printlnToFile(fout,ii,R2);
			
		}
		
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
	
	private double getR2(double[][] dataArray){
		
		double cnt = 0;		
		double cntR2 = 0;
		
		for (int jj = 0 ; jj < dataArray[0].length ; jj++){
			cnt+=dataArray[2][jj];
			cntR2+=dataArray[1][jj]*dataArray[1][jj]*dataArray[2][jj];
		}
		
		return cntR2/cnt;
	}

}

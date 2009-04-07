package chris.tests;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Random;
import chris.util.SortUtil;
import scikit.jobs.Control;
import scikit.jobs.Simulation;
import scikit.jobs.params.DirectoryValue;

public class sortTests extends Simulation{

	private String outdir;
	private double rlist[], slist[];
	private int idlist[];
	private static int N = 20;
	private Random rand;
	
	public static void main(String[] args) {
		new Control(new sortTests(), "Parameters");
	}
	
	public void animate() {
		
	}

	public void clear() {
		
	}

	public void load(Control c) {
		params.add("Data Directory",new DirectoryValue("/Users/cserino/Desktop/"));
		params.add("Seed",(int) 1);

	}

	public void run() {
		outdir = params.sget("Data Directory");
		rand   = new Random(params.iget("Seed"));
		
		rlist  = new double[N];
		slist  = new double[N];
		idlist = new int[N];
		
		for (int jj = 0 ; jj < N ; jj++){
			rlist[jj] = 5*rand.nextDouble();
			slist[jj] = rlist[jj];
		}
		
		idlist = SortUtil.S2LindexSort(slist);
		slist  = SortUtil.S2Lsort(slist);
	
		String fout = outdir + File.separator+"Sort.txt";
		
		try{
			File file = new File(fout);
			PrintWriter pw = new PrintWriter(new FileWriter(file, true), true);

			pw.print("Initial List");
			pw.print("\t");
			pw.print("Inverted Sorted");
			pw.print("\t");
			pw.print("Inverted Sorted via Indices");
			pw.print("\t");
			pw.print("Indices");
			pw.print("\t");
			pw.print("Difference");
			pw.println();
			pw.println();
			
			for (int jj = 0 ; jj < N ; jj++){		
				pw.print(rlist[jj]);
				pw.print("\t");
				pw.print(slist[jj]);
				pw.print("\t");
				pw.print(rlist[idlist[jj]]);
				pw.print("\t");
				pw.print(idlist[jj]);
				pw.print("\t");
				pw.print(slist[jj]-rlist[idlist[jj]]);
				pw.println();
				pw.println();
			}
			
			pw.close();
		}
		catch (IOException ex){
			ex.printStackTrace();
		}

		
		
	}
	
	
}

package chris.queue.Apps;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.text.DecimalFormat;
import java.util.LinkedList;

import chris.queue.message;
import chris.util.PrintUtil;
import chris.util.Random;
import scikit.jobs.Control;
import scikit.jobs.Job;
import scikit.jobs.Simulation;
import scikit.jobs.params.DirectoryValue;
import scikit.jobs.params.DoubleValue;

public class singleRouterApp extends Simulation{

	private LinkedList<message> buffer;
	private double lambda;
	private int L, Nmsg, tmax, dl, data[], Ns[];
	private Random rand;
	private String outdir, bname;
	
	public static void main(String[] args) {
		new Control(new singleRouterApp(), "Single Router");
	}
	
	public void animate() {
		
		params.set("messages",Nmsg);
		return;
	}

	public void clear() {
		params.set("messages",0);
		return;
	}

	public void load(Control c) {
		params.add("Data Directory",new DirectoryValue("/Users/cserino/Desktop/"));
		params.add("Data File", "default");
		params.add("l",5);
		params.add("\u03BB",new DoubleValue(0.05,0,1));
		params.add("seed",0);
		params.add("t_max",(int)(1e6));
		params.add("messages");
		params.set("messages",0);
		params.add("t");
		params.set("t",0);
	}

	public void run() {
		
		int idx, s, cycle;
		
		cycle = 0;
	
		L      = params.iget("l");
		rand   = new Random(params.iget("seed"));
		lambda = params.fget("\u03BB");
		outdir = params.sget("Data Directory");
		bname  = params.sget("Data File");
		tmax   = params.iget("t_max");
		Nmsg   = 0;
		dl     = L*tmax;
		data   = new int[dl];
		Ns     = new int[L];
		buffer = new LinkedList<message>();


		for(int jj = 0 ; jj < L ; jj++){
			data[jj] = 0;
			Ns[jj]   = 0;
		}

		PrintUtil.printlnToFile(outdir+File.separator+"Params_"+bname+".log",params.toString());

		while(true){

			for (int jj = 1 ; jj < tmax ; jj++){

				// select a message at random and "pass" it 
				if(buffer.size() > 0){
					idx = rand.nextInt(buffer.size());
					// get the hop number, pop it, and dissipate it if appropriate
					s = buffer.get(idx).getHops();
					if(s == L - 1){
						// dissipate
						buffer.remove(idx);
						Nmsg--;
						Ns[s]--;
					}
					else{
						// 
						buffer.get(idx).hopped();
						Ns[s]--;
						Ns[s+1]++;
					}
				}

				// try and generate a message
				if (rand.nextDouble() < lambda){
					buffer.add(new message(jj,0));
					Nmsg++;
					Ns[0]++;
				}

				// save distribution
				for(int kk = 0 ; kk < L ; kk++){
					data[kk + jj*L] = Ns[kk];  
				}

				if (jj % 1000 == 0){
					Job.animate();
				}

			}
			saveData(cycle);
			cycle++;
		}	
	}
	
	public void saveData(int ccl){
		
		DecimalFormat fmtS = new DecimalFormat("00");
		DecimalFormat fmtL = new DecimalFormat("0000000");

		try{
			File file = new File(outdir+File.separator+bname+fmtL.format(ccl)+".txt");
			PrintWriter pw = new PrintWriter(new FileWriter(file, true), true);
			pw.print("L = "+fmtS.format(L));
			for (int jj = 0 ; jj < tmax ; jj++){
				pw.print(jj);
				pw.print("\t");
				pw.println(data[jj]);
			}			
			pw.close();
		}
		catch (IOException ex){
			ex.printStackTrace();
		}
		return;
	}
	
}

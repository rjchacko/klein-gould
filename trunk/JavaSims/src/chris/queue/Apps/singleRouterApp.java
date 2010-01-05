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
import scikit.dataset.DatasetBuffer;
import scikit.dataset.Histogram;
import scikit.jobs.Control;
import scikit.jobs.Job;
import scikit.jobs.Simulation;
import scikit.jobs.params.DirectoryValue;
import scikit.jobs.params.DoubleValue;

public class singleRouterApp extends Simulation{

	private LinkedList<message> buffer;
	private double lambda;
	private int L, Nmsg, dl, data[], Ns[], now;
	private Random rand;
	private String outdir, bname;
	
	public static void main(String[] args) {
		new Control(new singleRouterApp(), "Single Router");
	}
	
	public void animate() {
		
		params.set("messages",Nmsg);
		params.set("t",now);
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
		params.add("t_ss",(int)(1e6));
		params.add("t_max",(int)(1e6));
		params.add("messages");
		params.set("messages",0);
		params.add("t");
		params.set("t",0);
		params.add("cycle");
		params.set("cycle",0);
	}

	public void run() {
		
		int idx, s, cycle, tss, tmax;
		Histogram h;
		
		cycle = 0;
	
		L      = params.iget("l");
		lambda = params.fget("\u03BB");
		outdir = params.sget("Data Directory");
		bname  = params.sget("Data File");
		tmax   = params.iget("t_max");
		tss    = params.iget("t_ss");
		dl     = L*tss;



		PrintUtil.printlnToFile(outdir+File.separator+"Params_"+bname+".log",params.toString());

		while(cycle < 1e5){

			buffer = new LinkedList<message>();
			rand   = new Random(params.iget("seed"));
			data   = new int[dl];
			Ns     = new int[L];
			Nmsg   = 0;

			// etc.
			
			
			// first simulate the approach to the steady state			
			for (int jj = 0 ; jj < tss ; jj++){

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
					data[kk+jj*L] = Ns[kk];  
				}

				if (jj % 1000 == 0){
					now = jj - tss;
					Job.animate();
				}

			}
			saveData(cycle, tss);
			
			// now simulate the (assumed) steady state
			h = new Histogram(1);
			for (int jj = 0 ; jj < tmax ; jj++){

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

				h.accum(Nmsg);
				
				if (jj % 1000 == 0){
					now = jj;
					Job.animate();
				}

			}
			saveHist(cycle, h);
			cycle++;
			params.set("cycle",cycle);;
			params.set("seed",params.iget("seed")+1);
		}	
	}
	
	public void saveData(int ccl, int tub){
		
		DecimalFormat fmtS = new DecimalFormat("00");
		DecimalFormat fmtL = new DecimalFormat("0000000");

		try{
			File file = new File(outdir+File.separator+bname+fmtL.format(ccl)+".txt");
			PrintWriter pw = new PrintWriter(new FileWriter(file, true), true);
			pw.println("L = "+fmtS.format(L));
			for (int jj = 0 ; jj < L*tub ; jj++){
				pw.println(data[jj]);
			}			
			pw.close();
		}
		catch (IOException ex){
			ex.printStackTrace();
		}
		return;
	}
	
	public void saveHist(int ccl, Histogram h){

		DecimalFormat fmtL = new DecimalFormat("0000000");	
		DatasetBuffer hh   = h.copyData();

		try{
			File file = new File(outdir+File.separator+bname+"_hist"+fmtL.format(ccl)+".txt");
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
		return;
	}

}

package chris.queue.Apps;

import java.io.File;
import java.text.DecimalFormat;
import java.util.LinkedList;

import scikit.dataset.Histogram;
import scikit.jobs.Control;
import scikit.jobs.Job;
import scikit.jobs.Simulation;
import scikit.jobs.params.DirectoryValue;
import scikit.jobs.params.DoubleValue;
import chris.queue.message;
import chris.util.PrintUtil;
import chris.util.Random;

public class singleRouterSusceptibilityApp extends Simulation{

	private LinkedList<message> buffer;
	private double lambda;
	private int L, N, now, cycle;
	private Random rand;
	private String outdir, bname;
	DecimalFormat fmt  = new DecimalFormat("000");

	
	public static void main(String[] args) {
		new Control(new singleRouterSusceptibilityApp(), "Single Router");
	}
	
	public void animate() {
		
		params.set("messages",N);
		params.set("t",now);
		return;
	}

	public void clear() {
		
		params.set("t", 0);
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
		params.add("cycle");
		params.set("cycle",0);
	}

	public void run() {
		
		int idx, s, tss, tmax;
		Histogram hN, td;
		double dlambda = 0.05;
		boolean ss;
		
		L      = params.iget("l");
		lambda = params.fget("\u03BB");
		outdir = params.sget("Data Directory");
		bname  = params.sget("Data File");
		tmax   = params.iget("t_max");
		tss    = 0;
		PrintUtil.printlnToFile(outdir+File.separator+"Params_"+bname+".log",params.toString());

		while(lambda < 0.2){
			cycle = 0;
			while(cycle++ < 20){

				buffer = new LinkedList<message>();
				rand   = new Random(params.iget("seed"));
				N      = 0;
				ss     = false;

				// first simulate the approach to the steady state			
				while(!ss){
					tss++;
					// select a message at random and "pass" it 
					if(buffer.size() > 0){
						idx = rand.nextInt(buffer.size());
						// get the hop number, pop it, and dissipate it if appropriate
						s = buffer.get(idx).getHops();
						if(s == L - 1){
							// dissipate
							buffer.remove(idx);
							N--;
						}
						else{
							// 
							buffer.get(idx).hop();
						}
					}

					// try and generate a message
					if (rand.nextDouble() < lambda){
						buffer.add(new message(tss));
						N++;
					}
					
					if(tss % 100000 == 0){
						now = -tss;
						Job.animate();
					}
					
					ss = (N > 1/(1-lambda*L));
				}

				// now simulate for 1e6 more time steps
				for(int jj = 0 ; jj < 1e6 ; jj++){
					tss++;
					// select a message at random and "pass" it 
					if(buffer.size() > 0){
						idx = rand.nextInt(buffer.size());
						// get the hop number, pop it, and dissipate it if appropriate
						s = buffer.get(idx).getHops();
						if(s == L - 1){
							// dissipate
							buffer.remove(idx);
							N--;
						}
						else{
							// 
							buffer.get(idx).hop();
						}
					}

					// try and generate a message
					if (rand.nextDouble() < lambda){
						buffer.add(new message(tss));
						N++;
					}
					
					if(tss % 100000 == 0){
						now = -tss;
						Job.animate();
					}					
				}
				
				// now simulate the (assumed) steady state
				hN = new Histogram(1);
				td = new Histogram(1);

				for (int jj = 0 ; jj < tmax ; jj++){
					// select a message at random and "pass" it 
					if(buffer.size() > 0){
						idx = rand.nextInt(buffer.size());
						// get the hop number, pop it, and dissipate it if appropriate
						s = buffer.get(idx).getHops();
						if(s == L - 1){
							// dissipate
							td.accum(jj+tss-buffer.get(idx).getTcreate());
							buffer.remove(idx);
							N--;
						}
						else{
							// hop
							buffer.get(idx).hop();
						}
					}

					// try and generate a message
					if (rand.nextDouble() < lambda){
						buffer.add(new message(jj+tss));
						N++;
					}

					hN.accum(N);
					if (jj % 100000 == 0){
						now = jj;
						Job.animate();
					}
				}
				saveHist(hN,td);
				params.set("seed",params.iget("seed")+1);
			}
		lambda += dlambda;
		}
	}	

	
	private void saveHist(Histogram h1, Histogram h2){
		
		PrintUtil.printHistToFile(outdir+File.separator+"n"+fmt.format(1000*lambda)+"_"+fmt.format(cycle)+".txt", h1);
		PrintUtil.printHistToFile(outdir+File.separator+"t"+fmt.format(1000*lambda)+"_"+fmt.format(cycle)+".txt", h2);
		return;
	}
}

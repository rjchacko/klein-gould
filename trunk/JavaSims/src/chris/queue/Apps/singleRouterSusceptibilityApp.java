package chris.queue.Apps;

import java.io.File;
import java.text.DecimalFormat;
import java.util.LinkedList;

import scikit.dataset.Histogram;
import scikit.jobs.Control;
import scikit.jobs.Job;
import scikit.jobs.Simulation;
import scikit.jobs.params.DirectoryValue;
import chris.queue.message;
import chris.util.PrintUtil;
import chris.util.Random;

public class singleRouterSusceptibilityApp extends Simulation{

	private LinkedList<message> buffer;
	private double lambda;
	private int L, N, now, cycle, p;
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
		//params.add("\u03BB",new DoubleValue(0.05,0,1));
		params.add("seed",0);
		params.add("t_sim",(int)(1e7));
		params.add("messages");
		params.set("messages",0);
		params.add("t");
		params.set("t",0);
		params.add("cycle");


	}

	public void run() {
		
		int idx, s, tss, tmax;
		Histogram hN, td;
		//double dlambda = 0.001;
		boolean ss;
		
		p      = 1;
		L      = params.iget("l");
		lambda = 0.1;
		outdir = params.sget("Data Directory");
		bname  = params.sget("Data File");
		tmax   = params.iget("t_sim");
		tss    = 0;
		PrintUtil.printlnToFile(outdir+File.separator+"Params_"+bname+".log",params.toString());

		while(p < 10){
			cycle = 0;
			System.out.println("Starting New Cycle");
			while(cycle++ < 10){
				System.out.println("Here 1");
				params.set("cycle",cycle);
				buffer = new LinkedList<message>();
				rand   = new Random(params.iget("seed"));
				N      = 0;
				ss     = false;
				System.out.println("Here 2");

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
					
					if(tss % 500000 == 0){
						now = (int)(100/(N*(1-lambda*L)));
						Job.animate();
					}
					
					ss = (N > 1/(1-lambda*L));
				}
				System.out.println("Here 3");

				// now simulate for 1e7 more time steps
				for(int jj = 0 ; jj < 1e7 ; jj++){
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
					
					if(jj % 500000 == 0){
						now = (int)(jj-1e7);
						Job.animate();
					}					
				}
				System.out.println("Here 4");

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
					if (jj % 500000 == 0){
						now = jj;
						Job.animate();
					}
				}
				saveHist(hN,td);
				params.set("seed",params.iget("seed")+1);
			}
			System.out.println("Next p");
		lambda = 2. - Math.pow(10,-p);
		p++;
		}
	}	

	
	private void saveHist(Histogram h1, Histogram h2){
		
//		PrintUtil.printHistToFile(outdir+File.separator+"n"+fmt.format(1000*lambda)+"_"+fmt.format(cycle)+".txt", h1);
//		PrintUtil.printHistToFile(outdir+File.separator+"t"+fmt.format(1000*lambda)+"_"+fmt.format(cycle)+".txt", h2);
		PrintUtil.printHistToFile(outdir+File.separator+"n"+fmt.format(p)+"_"+fmt.format(cycle)+".txt", h1);
		PrintUtil.printHistToFile(outdir+File.separator+"t"+fmt.format(p)+"_"+fmt.format(cycle)+".txt", h2);
		return;
	}
}

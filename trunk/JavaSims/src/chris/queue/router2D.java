package chris.queue;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.LinkedList;

import scikit.jobs.params.Parameters;
import chris.util.PrintUtil;
import chris.util.Random;

public class router2D {

	public static final int dl = 5000; 
	public static final int dc = 2;
	protected LinkedList<LinkedList<message>> buffer;
	protected double lambda, data[][], nbar[], omega;
	protected int N, LL, L, Nmsg, t;
	protected Random rand;
	protected String outdir, bname;
	
	public router2D(Parameters params){
		
		constructor_router2D(params);
		return;
	}
	
	public void constructor_router2D(Parameters params){
		
		t      = 0;
		LL     = params.iget("L");
		N      = LL*LL;
		L      = params.iget("l");
		lambda = params.fget("\u03BB");
		outdir = params.sget("Data Directory");
		bname  = params.sget("Data File");
		
		rand   = new Random(params.iget("seed"));
		buffer = new LinkedList<LinkedList<message>>(); // a list of the list's of buffers
		Nmsg   = 0;
		nbar   = new double[N];
		data   = new double[dc][dl];
		
		PrintUtil.printlnToFile(outdir+File.separator+"Params_"+bname+".log",params.toString());
		
		// create a list of buffers
		for (int jj = 0 ; jj < N ; jj++){
			buffer.add(new LinkedList<message>());
		}
		
		return;
	}
	
	
	public int step(boolean rec){
		
		return step(1, rec);
	}
	
	public int step(int ns, boolean rec){


		message[] tomove = new message[N];
		int idx, h;
		int tcount  = 0;
		double tbar = 0;
		
		for (int jj = 0 ; jj < ns ; jj++){
		
			t++;
			
			double tmp = 0;
			omega = 0;

			for (int kk = 0 ; kk < N ; kk++){
				// use loop to calculate metric (PART I)
				if(rec){
					nbar[kk] += buffer.get(kk).size();
					tmp      += nbar[kk];
				}

				// select a message to move at random from the buffer
				if(buffer.get(kk).isEmpty()) 
					continue;
				idx = rand.nextInt(buffer.get(kk).size());
				tomove[kk] = buffer.get(kk).get(idx);
				buffer.get(kk).remove(idx);
			}
			tmp /= N;
			for (int kk = 0 ; kk < N ; kk++){
				// use loop to calculate metric (PART II)
				if(rec)
					omega += (nbar[kk]-tmp)*(nbar[kk]-tmp);

				// generate new messages and add them to the buffers
				if (rand.nextDouble() < lambda){
					buffer.get(kk).add(new message(t,L,LL,kk,rand));
					Nmsg++;
				}

				// move messages selected in previous loop to next router
				if(tomove[kk] == null) 
					continue;
				h = tomove[kk].hop();
				if(h == L-1){	// dissipate message
					tbar += t-tomove[kk].getTcreate();
					tcount++;
					tomove[kk] = null; // clear from memory
					Nmsg--;
				}
				else{	// pass message to next router
					buffer.get(tomove[kk].getMove(h)).add(tomove[kk]);
				}
			}
		}

		// finish metric calculation
		if(rec){
			omega /= ((double)(N)*(double)(t)*(double)(t));
			takedata(omega, tbar/tcount);
		}
		
		return Nmsg;
	}
	
	public void takedata(double metric, double tb){

		if(t%dl == 0)
			writeData(t);
		
		data[0][t%dl] = 1/metric;
		data[1][t%dl] = tb;
		return;
	}
	
	public void writeData(int tt){
		
		int ub;
		int offset = (int)((tt-1)/dl);
		offset = offset*dl;

		ub = tt%dl;
		if(ub==0) ub = dl;
		
		try{
			File file = new File(outdir+File.separator+bname+".txt");
			PrintWriter pw = new PrintWriter(new FileWriter(file, true), true);
			for (int jj = 0 ; jj < ub ; jj++){
				pw.print(jj+offset);
				pw.print("\t");
				for (int kk = 0 ; kk < dc ; kk++){
					pw.print(data[kk][jj]);
					pw.print("\t");
				}		
				pw.println();
			}
			pw.close();
		}
		catch (IOException ex){
			ex.printStackTrace();
		}
		return;
	}
	
	public int[] getDensity(){
		
		int[] ret = new int[N];
		
		for (int jj = 0 ; jj < N ; jj++){
			ret[jj] = buffer.get(jj).size();
		}
		
		return ret;
	}
	
	public int getT(){
		
		return t;
	}
	
	public int getNmsg(){
		
		return Nmsg;
	}

}

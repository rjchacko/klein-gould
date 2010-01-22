package chris.queue;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.text.DecimalFormat;

import scikit.jobs.params.Parameters;


public class finiteRouter2D extends router2D{
	
	private int M, npp[], nad[], ccl;
	public static final DecimalFormat cmt = new DecimalFormat("000");

	
	public finiteRouter2D(Parameters params) {
		
		super(params);
		constructor_finiteRouter2D(params);
	}

	public void constructor_finiteRouter2D(Parameters params){
	
		M   = params.iget("M");
		nad = new int[dl*2];
		npp = new int[dl*N];
		ccl = params.iget("cycle");
		
		return;
	}
	
	public int step(int ns, boolean rec){

		message[] tomove = new message[N];
		int idx, h;
		int tcount  = 0;
		double tbar = 0;
		int na = 0;
		int nd = 0;
		
		for (int jj = 0 ; jj < ns ; jj++){
			
			t++;
			if(t%dl == 0)
				writePPdata(getT(), ccl);
				
				
			double tmp = 0;
			omega = 0;

			for (int kk = 0 ; kk < N ; kk++){
				// use loop to calculate metric (PART I)
				if(rec){
					nbar[kk] += buffer.get(kk).size();
					tmp      += nbar[kk];
				}
				
				// generate new messages and add them to the buffers
				if (rand.nextDouble() < lambda && buffer.get(kk).size() < M){	// message generated
					buffer.get(kk).add(new message(t,L,LL,kk,rand));
					na++;
					Nmsg++;
					
					// select messages to move
					if(buffer.get(kk).size() == 1)	
						continue;	// no message to move
					idx = rand.nextInt(buffer.get(kk).size()-1); // do not choose message just made
					if(buffer.get(buffer.get(kk).get(idx).getMove()).size() >= M )
						continue;	// destination is full
						
					tomove[kk] = buffer.get(kk).get(idx);
				}
				else{	// no message generated
					// select messages to move
					if(buffer.get(kk).isEmpty())	
						continue;	// no message to move
					idx = rand.nextInt(buffer.get(kk).size()); // choose from all messages
					if(buffer.get(buffer.get(kk).get(idx).getMove()).size() >= M )
						continue;	// destination is full
						
					tomove[kk] = buffer.get(kk).get(idx);
				}
			}
			for (int kk = 0 ; kk < N ; kk++){
				// use loop to calculate metric (PART II)
				if(rec)
					omega += (nbar[kk]-tmp)*(nbar[kk]-tmp);

				// move messages selected in previous loop to next router
				if(tomove[kk] == null || buffer.get(tomove[kk].getMove()).size() >= M) 
					continue; // no message to move OR destination has become full
			
				buffer.get(kk).remove(tomove[kk]);
				h = tomove[kk].hop();
				if(h == L-1){	// dissipate message
					tbar += t-tomove[kk].getTcreate();
					tcount++;
					tomove[kk] = null; // clear from memory
					nd++;
					Nmsg--;
				}
				else{	// pass message to next router
					buffer.get(tomove[kk].getMove(h)).add(tomove[kk]);
				}
				npp[(kk+N*t)%dl] = buffer.get(kk).size();
			}
		}
		
		// finish metric calculation
		if(rec){
			omega /= ((double)(N)*(double)(t)*(double)(t));
			takedata(omega, tbar/tcount);
		}
		nad[(2*t)%dl]   = na;
		nad[(1+2*t)%dl] = nd;
		
		return Nmsg;
	}
	
	public void writePPdata(int tt, int cycle){
		
		int ub;
		int offset = (int)((tt-1)/dl);
		offset = offset*dl;

		ub = tt%dl;
		if(ub==0) ub = dl;
		
		try{
			File file = new File(outdir+File.separator+bname+"_NofT"+cmt.format(1000*lambda)+"_"+cmt.format(cycle)+".txt");
			PrintWriter pw = new PrintWriter(new FileWriter(file, true), true);
			pw.println("N = " + N);
			for (int jj = 0 ; jj < ub ; jj++)
				for (int kk = 0 ; kk < N ; kk++)
					pw.println(npp[kk+jj*N]);
			pw.close();
		}
		catch (IOException ex){
			ex.printStackTrace();
		}
		try{
			File file = new File(outdir+File.separator+bname+"_Nad"+cmt.format(1000*lambda)+"_"+cmt.format(cycle)+".txt");
			PrintWriter pw = new PrintWriter(new FileWriter(file, true), true);
			for (int jj = 0 ; jj < ub ; jj++){
				pw.print(jj+offset);
				pw.print("\t");
				pw.print(nad[2*jj]);
				pw.print("\t");
				pw.println(nad[2*jj+1]);
			}
			pw.close();
		}
		catch (IOException ex){
			ex.printStackTrace();
		}
		
		
		return;
	}
}

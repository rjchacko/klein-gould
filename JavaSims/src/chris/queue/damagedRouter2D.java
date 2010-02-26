package chris.queue;

import chris.util.SortUtil;
import scikit.dataset.Histogram;
import scikit.jobs.params.Parameters;

public class damagedRouter2D extends router2D{

	protected double p;
	protected boolean[] dn;
	
	public damagedRouter2D(Parameters params){
		
		super(params);
		damagedRouter2D_contructor(params);
	}
	
	public void damagedRouter2D_contructor(Parameters params){
		
		int[] rs;
		p             = params.fget("p");		
		int Ndn       = (int)(Math.round(p*N));
		dn            = new boolean[N];
		double rseq[] = new double[N];
	
		for (int jj = 0 ; jj < N ; jj++){
			rseq[jj] = rand.nextDouble();
			dn[jj]   = false;   
		}
		rs = SortUtil.S2LindexSort(rseq);
		
		for(int jj = 0 ; jj < Ndn ; jj++){
			dn[rs[jj]] = true;
		}
		
		return;
	}
	
	
	public int step(int ns, boolean takeData, Histogram hh){


		message[] tomove = new message[N];
		int idx, h;		

		
		for (int jj = 0 ; jj < ns ; jj++){
			
			int na = 0;
			int nd = 0;
			int tcount  = 0;
			double tbar = 0;
			
			t++;
			if(takeData)
				tm++;
			
			double tmp = 0;
			omega = 0;

			for (int kk = 0 ; kk < N ; kk++){
				if(dn[kk])
					continue;
				// use loop to calculate metric (PART I)
				if(takeData){
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
				if(takeData)
					omega += (nbar[kk]-tmp)*(nbar[kk]-tmp);

				// generate new messages and add them to the buffers
				if (rand.nextDouble() < lambda){
					buffer.get(kk).add(new message(t,L,LL,kk,rand));
					Nmsg++;
					na++;
				}

				// move messages selected in previous loop to next router
				if(tomove[kk] == null) 
					continue;
				h = tomove[kk].hop();
				if(h == L-1 || dn[tomove[kk].getMove(h)]){	// dissipate message
					tbar += t-tomove[kk].getTcreate();
					hh.accum(t-tomove[kk].getTcreate());
					tcount++;
					tomove[kk] = null; // clear from memory
					Nmsg--;
					nd++;
				}
				else{	// pass message to next router
					buffer.get(tomove[kk].getMove(h)).add(tomove[kk]);
				}
			}
			// finish metric calculation
			if(takeData){
				omega /= ((double)(N)*(double)(tm)*(double)(tm));
				takedata(omega, tbar/tcount, Nmsg, na, nd);
			}			
		}

		
		return Nmsg;
	}

}

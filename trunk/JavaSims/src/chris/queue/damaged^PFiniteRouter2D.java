package chris.queue;

import scikit.dataset.Histogram;
import scikit.jobs.params.Parameters;

public class damagedFiniteRouter2D extends damagedRouter2D{
	
	private int M;

	public damagedFiniteRouter2D(Parameters params){
		
		super(params);
		damagedFiniteRouter2D_contructor(params);
	}

	public void damagedFiniteRouter2D_contructor(Parameters params){
		
		M = params.iget("M");
		return;
	}
	
	public int step(int ns, boolean takeData, Histogram hh){
		/// add hist collection
		message[] tomove = new message[N];
		int idx, h;
		
		for (int jj = 0 ; jj < ns ; jj++){
			
			int tcount  = 0;
			double tbar = 0;
			int na = 0;
			int nd = 0;
			
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
			tmp /= N;
			for (int kk = 0 ; kk < N ; kk++){
				if(dn[kk])
					continue;
				// use loop to calculate metric (PART II)
				if(takeData)
					omega += (nbar[kk]-tmp)*(nbar[kk]-tmp);

				// move messages selected in previous loop to next router
				if(tomove[kk] == null || buffer.get(tomove[kk].getMove()).size() >= M) 
					continue; // no message to move OR destination has become full
			
				buffer.get(kk).remove(tomove[kk]);
				h = tomove[kk].hop();
				if(h == L-1 || dn[tomove[kk].getMove(h)]){	// dissipate message
					tbar += t-tomove[kk].getTcreate();
					tcount++;
					tomove[kk] = null; // clear from memory
					nd++;
					Nmsg--;
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
	
	public void setM(int M){
		
		this.M = M;
	}
}

package chris.queue;

import scikit.jobs.params.Parameters;


public class finiteRouter2D extends router2D{
	
	private int M;
	
	public finiteRouter2D(Parameters params) {
		
		super(params);
		constructor_finiteRouter2D(params);
	}

	public void constructor_finiteRouter2D(Parameters params){
	
		M = params.iget("M");
		return;
	}
	
	public void step(int ns, boolean rec){

		message[] tomove = new message[N];
		int idx, h;
		int tcount  = 0;
		double tbar = 0;
		
		
		for (int jj = 0 ; jj < ns ; jj++){
			
			if(rec)
				t++;
			
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
		return;
	}
}

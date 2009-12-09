package chris.queue;

import java.util.LinkedList;

import scikit.jobs.params.Parameters;
import chris.util.Random;

public class router1D {

	private LinkedList<LinkedList<message>> buffer;
	private double lambda;
	private int N, L, Nmsg;
	private int t;
	private Random rand;
	
	public router1D(Parameters params){
		
		constructor_router1D(params);
		return;
	}
	
	public void constructor_router1D(Parameters params){
		
		t      = 0;
		N      = params.iget("N");
		L      = params.iget("l");
		lambda = params.fget("\u03BB");
		rand   = new Random(params.iget("seed"));
		buffer = new LinkedList<LinkedList<message>>();
		Nmsg   = 0;
		
		// create a list of buffers
		for (int jj = 0 ; jj < N ; jj++){
			buffer.add(new LinkedList<message>());
		}
		
		return;
	}
	
	public void step(){
		
		step(1);
		return;
	}
	
	public void step(int ns){
		
		for (int jj = 0 ; jj < ns ; jj++){

			message[] tomove = new message[N];
			int idx, nidx;
			double r;

			t++;
			
			// Visit each router, select a message at random
			// from its buffer, and move it to the next router's buffer
			for (int kk = 0 ; kk < N ; kk++){
				if(buffer.get(kk).isEmpty()) 
					continue;
				idx = rand.nextInt(buffer.get(kk).size());
				tomove[kk] = buffer.get(kk).get(idx);
				tomove[kk].hopped();
				buffer.get(kk).remove(idx);
			}
			for (int kk = 0 ; kk < N ; kk++){
				if(tomove[kk] == null) 
					continue;
				//System.out.println(tomove[kk].getHops());
				if(tomove[kk].getHops() == L){
					// dissipate message
					tomove[kk] = null; // clear from memory
					Nmsg--;
				}
				else{
					// pass message to next router
					nidx = kk + tomove[kk].getMove(0);
					if(nidx == -1 || nidx == N) 
						nidx = ((N-nidx)/(N+1))*(N-1); // maps -1 --> N-1 , N --> 0 
					buffer.get(nidx).add(tomove[kk]);
				}
			}
			
			// generate new messages and add them to the buffers
			for (int kk = 0 ; kk < N ; kk++){
				r = rand.nextDouble();
				if (r < lambda){
					r = r < lambda/2 ? -1 : 1;
					buffer.get(kk).add(new message(t,(int)(r)));
					Nmsg++;
				}
			}
			
//			/*
//			 *  DEBUGGING
//			 */
//			if (t==1){
//	
//				buffer.get(0).add(new message(t,1));
//				buffer.get(4).add(new message(t,-1));
//				
//				Nmsg++;
//				Nmsg++;
//			}
//			if (t==3){
//				
//				buffer.get(2).add(new message(t,1));
//				
//				Nmsg++;
//				Nmsg++;
//			}
//			
//			/*
//			 * DEBUGGING
//			 * 
//			 */
			
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

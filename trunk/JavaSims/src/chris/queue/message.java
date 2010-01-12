package chris.queue;

import chris.util.Random;
import chris.util.SortUtil;


public class message {
	
	private int tc, dir[], hops;
	
	public message(int tc){
		
		constructor_message(tc);
		return;
	}
	
	public message(int tc,  int L, int N, int x, int dir){
		
		constructor_message(tc, L, N, x, dir);
		return;
	}
	
	public message(int tc, int L, int LL, int r, Random rand){
		
		constructor_message(tc, L, LL, r, rand);
		return;
	}
	
	// Zero-dimensional (i.e. single buffer) geometry constructor
	public void constructor_message(int tc){
		
		this.tc = tc;
		dir     = null;
		hops    = 0;
	}
		
	// Ring geometry constructor
	public void constructor_message(int tc, int L, int N, int x, int dir){
		

		hops        = 0;
		this.dir    = new int[L];
		this.tc     = tc;
		
		if(dir < 0){
			// find end point and walk backward
			int ep = x-L;
			for (int jj = 0 ; jj < L ; jj++){
				this.dir[L-1-jj] = (ep+N+jj)%N;
			}
		
		}
		else{
			// walk forward 
			for (int jj = 0 ; jj < L ; jj++){
				this.dir[jj] = (x+jj+1)%N;
			}
		}
	}
	
	// Torus nearest neighbor geometry constructor
	public void constructor_message(int tc, int L, int LL, int r, Random rand){
		
		int dest, dx, dy, rx, ry;

		hops          = 0;
		this.tc       = tc;
		dir           = new int[L];
		boolean[] tmp = new boolean[L];
		rx            = r%LL;
		ry            = (int)(r/LL);
		
		// generate destination and route for message
		dest = rand.nextInt(4*L);
		dx   = dest > L ? Math.abs(dest-3*L)-L : dest;
		dy   = Math.abs(dest-2*L)-L;
		for(int jj = 0 ; jj < Math.abs(dx) ; jj++){
			tmp[jj] = true;
		}
		for(int jj = 0 ; jj < Math.abs(dy) ; jj++){
			tmp[jj+Math.abs(dx)] = false;
		}
		dx = (int)Math.signum(dx);
		dy = (int)Math.signum(dy);
		// Shuffle the order of the steps
		tmp = SortUtil.shuffleArray(tmp, rand);
		for (int jj = 0 ; jj < L ; jj++){
			if(tmp[jj]){	// step left / right
				rx      = (rx+dx+LL)%LL;
				dir[jj] = rx + ry*LL; 					
			}
			else{	// step up / down
				ry      = (ry+dy+LL)%LL;
				dir[jj] = rx + ry*LL; 
			}
		}
		return;
	}
	
	public int getMove(int step){
		
		return dir[step];
	}
	
	public int getMove(){
		
		return dir[hops];
	}
	
	public int getTcreate(){
		
		return tc;
	}
	
	public int getHops(){
		
		return hops;
	}
	
	public int hop(){
		
		return hops++;
	}
}

package chris.queue;


public class message {
	
	private int tc, dir[], dim, hops;
	
	public message(int tc, int dir){
		
		constructor_message(tc, dir);
		return;
	}
	
	public void constructor_message(int tc, int dir){
		
		hops        = 0;
		dim         = 1;
		this.dir    = new int[dim];
		this.tc     = tc;
		this.dir[0] = dir;
		
		return;
	}
	
	public int getMove(int step){
		
		return dir[step];
	}
	
	public int getTcreate(){
		
		return tc;
	}
	
	public int getHops(){
		
		return hops;
	}
	
	public void hopped(){
		
		hops++;
	}
	
}

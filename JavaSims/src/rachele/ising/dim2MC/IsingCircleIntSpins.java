package rachele.ising.dim2MC;

public class IsingCircleIntSpins {
	public int N, L, R;
	public int [] spin;
	int netSum;	

	public IsingCircleIntSpins(int L, int R) {
		this.L = L;
		this.R = R;
		N = L*L;
		netSum = N;
		spin = new int[L*L];
		for (int i = 0; i < L*L; i ++){
			spin[i] = 1;
		}
	}
	
	public int sumAll() {
		return netSum;
	}
	
	public int sumInRange(int x, int y) {
		int sum = 0;
		for (int xadd = -R; xadd <= R; xadd ++){
			int xx = x + xadd;
			for(int yadd = -R; yadd <= R; yadd ++){
				int yy = y + yadd;
				double distance = Math.sqrt(yadd*yadd+xadd*xadd);
				if (distance <= R){
					int newPos = ((yy + L)%L)*L + (xx + L)%L;
					sum += spin[newPos];
				}
			}
		}
		return sum;
	}
	
	public int get(int x, int y) {
		int i = y*L+x;
		return spin[i];
	}
	
	public int findZ(){
		int sum = 0;
		for (int xadd = -R; xadd <= R; xadd ++){
			for(int yadd = -R; yadd <= R; yadd ++){
				double distance = Math.sqrt(yadd*yadd+xadd*xadd);
				if (distance <= R){
					sum += 1;
				}
			}
		}
		return sum;
	}
	
	public int get(int i) {
		return spin[i];
	}	
	
	public void set(int x, int y, int s) {
		assert (s == 1 || s == -1);
		if (get(x, y) != s)
			flip(x, y);
	}
	
	public void flip(int x, int y) {
		int dm = -2*get(x,y); // change in magnetization after flipping spin i
		int i = L*y+x;
		spin[i] *=-1;
		netSum += dm;
	}
	
}

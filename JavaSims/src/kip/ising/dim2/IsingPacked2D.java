package kip.ising.dim2;

//import kip.util.Random;

public class IsingPacked2D {
/*
	public int spin[];
	private int packed[];
	
	public int L1, L2;
	public int N;
	public double T;
	public Random random = new Random();
	public double time;
	public boolean openBoundary = false;
	
	
	public IsingPacked2D(int seed, int _L1, int _L2, double _T) {
		if (L1 % 32 != 0) {
			throw new IllegalArgumentException("System width must be multiple of 32");
		}

		random.setSeed(seed);
		L1 = _L1;
		L2 = _L2;
		T = _T;
		N = L1*L2;
		time = 0;
		
		spin = new int[N];
		packed = new int[N/32];
		randomize();
	}
	
	public void randomize() {
		for (int i = 0; i < N/32; i++)
			packed[i] = random.nextInt();
		unpack();
	}
	
	public void unpack() {
		for (int i = 0; i < N/32; i++) {
			int x = (i*32) % L1;
			int y = (i*32) / L1;
			
			for (int j = 0; j < 32; j++)
				spin[y*L1+(x+j)] = 2 * ((packed[i] >>> j) & 1) - 1;
		}
	}
	
	private int neighborSum(int i) {
		int x = i % L1;
		int y = i / L1;
		int acc = 0;
		
		if (openBoundary) {
			if (y < L2-1)
				acc += spin[i+L1];
			if (y > 0)
				acc += spin[i-L1];
			if (x < L1-1)
				acc += spin[i+1];
			if (x > 0)
				acc += spin[i-1];
		}
		else {
			int yp = (y+1)%L2;
			int ym = (y-1+L2)%L2;
			int xp = (x+1)%L1;
			int xm = (x-1+L1)%L1;
			acc += spin[yp*L1+x];
			acc += spin[ym*L1+x];
			acc += spin[y*L1+xp];
			acc += spin[y*L1+xm];
		}
		
		return acc;
	}
	
	private void singleStep() {
		int i = random.nextInt(N);
		double dE = 2*spin[i]*neighborSum(i);
		if (dE <= 0 || (T != 0 && random.nextDouble() < Math.exp(-dE/T))) {
			spin[i] = -spin[i];
		}
	}
	
	
	
	public void getNeighbors(int i, int[] ret) {
		int x = i % (L1/32) + L1/32;
		int y = i / (L1/32) + L2;
		
		int xp = (x+1)%(L1/32);
		int xm = (x-1)%(L1/32);
		int yp = (y+1)%(L2);
		int ym = (y-1)%(L2);
		
		ret[0] = y*(L1/32) + xp;
		ret[1] = y*(L1/32) + xm;
		ret[2] = yp*(L1/32) + x;
		ret[3] = ym*(L1/32) + x;
	}
	
	private void updateEven() {
		int n[] = new int[4];
		
		int a1=0, a2=0;
		for (int i = 0; i < N/32; i++) {
			getNeighbors(i, n);
			
			int EVEN_BITS = 0x55555555; // = 0b01010101 ...
			
		}
	}
	
	
	public void step(double mcs) {
		long imcs = Math.round(mcs);
		for (int i = 0; i < imcs; i++) {
			update(0);
			update(1);
		}
		time += mcs;
		unpack();
	}
	*/
}

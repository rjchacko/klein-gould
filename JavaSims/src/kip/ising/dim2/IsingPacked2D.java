package kip.ising.dim2;

//import kip.util.Random;

public class IsingPacked2D {
/*	public int spin[];
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
	
	int indexSpin(int spins, int b) {
		return 2*((spins>>>(31-b))&1)-1;
	}
	
	public void unpack() {
		for (int i = 0; i < N/32; i++) {
			int x = (i*32) % L1;
			int y = (i*32) / L1;
			
			for (int b = 0; b < 32; b++)
				spin[y*L1+(x+b)] = indexSpin(packed[i], b);
		}
	}
	
	// unused
//	private int neighborSum(int i) {
//		int x = i % L1;
//		int y = i / L1;
//		int acc = 0;
//		
//		if (openBoundary) {
//			if (y < L2-1)
//				acc += spin[i+L1];
//			if (y > 0)
//				acc += spin[i-L1];
//			if (x < L1-1)
//				acc += spin[i+1];
//			if (x > 0)
//				acc += spin[i-1];
//		}
//		else {
//			int yp = (y+1)%L2;
//			int ym = (y-1+L2)%L2;
//			int xp = (x+1)%L1;
//			int xm = (x-1+L1)%L1;
//			acc += spin[yp*L1+x];
//			acc += spin[ym*L1+x];
//			acc += spin[y*L1+xp];
//			acc += spin[y*L1+xm];
//		}
//		
//		return acc;
//	}
	
	// unused
//	private void singleStep() {
//		int i = random.nextInt(N);
//		double dE = 2*spin[i]*neighborSum(i);
//		if (dE <= 0 || (T != 0 && random.nextDouble() < Math.exp(-dE/T))) {
//			spin[i] = -spin[i];
//		}
//	}
	
	
	
	public void getNeighbors(int i, int[] ret) {
		int x = i % (L1/32) + L1/32;
		int y = i / (L1/32) + L2;
		
		int xp = (x+1)%(L1/32);
		int xm = (x-1)%(L1/32);
		int yp = (y+1)%(L2);
		int ym = (y-1)%(L2);
		
		ret[0] = y*(L1/32) + xp; // right
		ret[1] = y*(L1/32) + xm; // left
		ret[2] = yp*(L1/32) + x; // up
		ret[3] = ym*(L1/32) + x; // down
	}
	
	private void updateEven() {
		int n[] = new int[4];
		
		int a1=0, a2=0;
		for (int i = 0; i < N/32; i++) {
			getNeighbors(i, n);
			int here = packed[i];
			int right = packed[n[0]];
			int left = packed[n[1]];
			int up = packed[n[2]];
			int down = packed[n[3]];
			
			int EVEN_BITS = 0x55555555; // = 0b0101...0101
			int ODD_BITS = EVEN_BITS << 1; // = 0b1010...1010
			int vert = EVEN_BITS&up + EVEN_BITS&down;
			int horiz = ((ODD_BITS&here)>>>1) + ((ODD_BITS&here)<<1) + (right>>>31); 
			
			
			for (int j = 1; j < 32; j += 2) {
				// shift: 30, 2, 0
				int shift = (31-j);
//				int neighborSum = (vert >>> shift) & 0x3
//				double dE = 2*indexSpin(packed[i], j)*neighborSum;
//			if (dE <= 0 || (T != 0 && random.nextDouble() < Math.exp(-dE/T))) {
//				spin[i] = -spin[i];
//			}
			}
			
		}
	}
	
	
	public void step(double mcs) {
		long imcs = Math.round(mcs);
		for (int i = 0; i < imcs; i++) {
//			update(0);
//			update(1);
		}
		time += mcs;
		unpack();
	}*/
}

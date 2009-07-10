package chris.TFB;

import java.util.Random;

import scikit.jobs.params.Parameters;

public class TFB {

	private boolean fibre[];
	private int Nin, N;
	private double stress, K, D, beta, phi, ebar[], Dp; //phi = 1 ==> all intact
	private Random rand;

	public TFB(Parameters params){
		
		TFB_contructor(params);	
	}
	
	public void TFB_contructor(Parameters params){
		int seed;
		
		seed   = params.iget("seed");
		N      = params.iget("N");
		stress = params.fget("stress / fibre"); // Proportional to m_{load}*g / N
		K      = params.fget("K");
		D      = params.fget("D");
		beta   = 1./params.fget("T");
		Dp     = stress*stress/(2*K*N);
		
		Nin    = N;
		phi    = 1.;
		rand   = new Random(seed);
		fibre  = new boolean[N];
		ebar   = new double[N];
		
		for(int jj = 0 ; jj < N ; jj++){
			fibre[jj] = true;
		}
		
	}

	public void nextBundle(){
		
		int s;
		
		for (int jj = 0 ; jj < N ; jj++){
			s = rand.nextInt(N);
			if(flip(s)){
				fibre[s] = !(fibre[s]);
				Nin += bool2sgn(fibre[s]);
			}	
			phi = (double)(Nin) / (double)(N);
//			// DEBUG
//			PrintUtil.printlnToFile("/Users/cserino/Desktop/test/debug.txt",phi);
		}
		
	}
	
	private int bool2sgn(boolean b){
		
		return b ? 1 : -1;
	}
	
	private boolean flip(int sj){
		
		if(Nin == 1 && fibre[sj]) return false; // dh --> +infinty
		double dh = (Dp/phi - D)*bool2sgn(fibre[sj]);
		if (dh <  0) return true;
		return (rand.nextDouble() < Math.exp(-beta*dh));
	}

	public double getOmegaInv(double mct){
		
		double Omega  = 0;
		double tmpbar = 0;
		
		for(int jj = 0 ; jj < N ; jj++){
			ebar[jj] += getHj(fibre[jj]);
			tmpbar   += ebar[jj];
		}
		tmpbar = tmpbar / N;
		for(int jj = 0 ; jj < N ; jj++){
			Omega += (ebar[jj]-tmpbar)*(ebar[jj]-tmpbar);
		}
		
		return (mct*mct*N)/Omega;
	}
	
	private double getHj(boolean sj){
		
		double strain = stress/(K*N*phi);
		return sj ? (0.5*K*strain*strain - D) : 0;
	}
	
	public double getPhi(){
		
		return phi;
	}
	
	public double getE(){
		
		double strain = stress/(K*N*phi);
		return (0.5*K*strain*strain-D)*phi*N;
	}
	
}


package chris.TFB;

import java.util.Random;

import scikit.jobs.params.Parameters;

public class TFB {

	private boolean fibre[];
	private int Nin, N;
	private double stress, K, D, beta, phi, ebar[], fbar[], sbar[], Dp, Omega, Omega2, Omega3; //phi = Nalive / N
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
		fbar   = new double[N];
		sbar   = new double[N];
		
		for(int jj = 0 ; jj < N ; jj++){
			fibre[jj] = true;
		}
		
	}

	public void nextBundle(){
		
		int s;
		
		//for (int jj = 0 ; jj < N ; jj++){
			s = rand.nextInt(N);
			if(flip(s)){
				fibre[s] = !(fibre[s]);
				Nin += bool2sgn(fibre[s]);
			}	
			phi = (double)(Nin) / (double)(N);
		//}
		
	}
	
	private int bool2sgn(boolean b){
		
		return b ? 1 : -1;
	}
	
	private boolean flip(int sj){
		
		if(Nin == 1 && fibre[sj]) return false; // dh --> +infinity
		double dh = (Dp/phi - D)*bool2sgn(fibre[sj]);
		if (dh <  0) return true;
		return (rand.nextDouble() < Math.exp(-beta*dh));
	}

	public void calcOmega(double mct){
		
		Omega  = 0;
		Omega2 = 0;
		Omega3 = 0;
		double tmpbar  = 0;
		double tmpbar2 = 0;
		double tmpbar3 = 0;
		
		for(int jj = 0 ; jj < N ; jj++){
			ebar[jj] += getHj(fibre[jj]);
			tmpbar   += ebar[jj];
			fbar[jj] += (bool2sgn(fibre[jj])+1)/(double)(2);
			tmpbar2  += fbar[jj];
			sbar[jj] += ((bool2sgn(fibre[jj])+1)/(double)(2))/phi;
			tmpbar3  += sbar[jj];
//			PrintUtil.printlnToFile("/Users/cserino/desktop/db.txt", getHj(fibre[jj]), (bool2sgn(fibre[jj])+1)/(double)(2));
			
		}
		tmpbar  = tmpbar / N;
		tmpbar2 = tmpbar2 / N;
		tmpbar3 = tmpbar3 / N;
		
//		PrintUtil.printlnToFile("/Users/cserino/desktop/db.txt", "---------------------------------");
		
		for(int jj = 0 ; jj < N ; jj++){
			Omega  += (ebar[jj]-tmpbar)*(ebar[jj]-tmpbar);
			Omega2 += (fbar[jj]-tmpbar2)*(fbar[jj]-tmpbar2);
			Omega3 += (sbar[jj]-tmpbar3)*(sbar[jj]-tmpbar3);
		}
//		
//		PrintUtil.printlnToFile("/Users/cserino/desktop/db.txt", Omega, Omega2);
//		PrintUtil.printlnToFile("/Users/cserino/desktop/db.txt", "___________________________");
//		
		Omega  = Omega/(mct*mct*N);
		Omega2 = Omega2/(mct*mct*N);
		Omega3 = Omega3/(mct*mct*N);
	}
	
	public double getOmega1inv(){
	
		return 1/Omega;
	}
	
	public double getOmega2inv(){
		
		return 1/Omega2;
	}
	
	public double getOmega3inv(){
		
		return 1/Omega3;
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


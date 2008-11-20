package rachele.ising.dim2;

import static java.lang.Math.min;
import rachele.ising.dim2.IsingLR.DynType;
import scikit.jobs.params.Parameters;
import kip.util.Random;

public class IsingArbitraryConnect{

	Random random = new Random();
	public DynType dynamics = DynType.GLAUBER;
	public int [] spin;
	public int L, R;
	int z, sumSpins;
	public double T,h;
	double J, time;
	int [][] neighborTable;

	public IsingArbitraryConnect(Parameters params) {
		L = Integer.highestOneBit(params.iget("L"));
		params.set("L", L);
		spin = new int [L*L];
		R = min(params.iget("R"), L/2-1);
		params.set("R", R);
		z = (2*R+1)*(2*R+1)-1;
		neighborTable = new int [L*L][z];
		random.setSeed(params.iget("Random seed"));
		setParameters(params);
	}
	
	public void setParameters(Parameters params) {
		String dyn = params.sget("Dynamics", "Ising Glauber");
		if (dyn.equals("Ising Glauber"))
			dynamics = DynType.GLAUBER;
		else if (dyn.equals("Ising Metropolis"))
			dynamics = DynType.METROPOLIS;
		else if (dyn.equals("Kawasaki Glauber"))
			dynamics = DynType.KAWA_GLAUBER;
		else if (dyn.equals("Kawasaki Metropolis"))
			dynamics = DynType.KAWA_METROPOLIS;
		

		T  = params.fget("T");
		J  = params.fget("J", 1);
		h  = params.fget("h", 0.0);

	}
	
	public void randomizeField(double m) {
		if (m == 1 || m == -1) {
			for (int i = 0; i < L*L; i++)
				spin[i]=(int)m;
			sumSpins = (int)m*L*L;
		}
		else {
			int sum = 0;
			for (int i = 0; i < L*L; i++) {
				// p(s = +-1) = (1 +- m) / 2
				int s = (random.nextDouble() < (1+m)/2) ? 1 : -1;
				spin[i]=s;
				sum += s;
			}
			sumSpins = sum;
		}
	}
	
	public void step() {
//		for(int i = 0; i < L*L; i++){
			//choose random spin
			int randomSpin =  (int)(random.nextDouble()*(double)(L*L));
			//calculate dE
			int nborSum = 0;
			for (int n = 0; n < z; n++)
				nborSum += spin[neighborTable[randomSpin][n]];
			double dE = 2*spin[randomSpin]*(h + J*(nborSum)/(z));
			//calculate prob of flipping for Glauber
			double prob = 1.0/(1.0+Math.exp(dE/T));
			//throw random number and choose flip
			if(random.nextDouble() < prob){
				spin[randomSpin] *= -1;
				sumSpins += 2*spin[randomSpin];
			}
//		}
		time += 1.0/(L*L);
	}
		
	
	public double magnetization() {
		return (double)sumSpins / (L*L);
	}
	
	public double time(){
		return time;
	}
	
	public int [] getField(int dx){
		return spin;
	}
	
	public void setSquareNeighbors(){	
		for(int i = 0; i < L*L; i++){
			int x = i%L;
			int y = i/L;
			int nborIndex = 0;
			for(int dx = -R; dx <= R; dx++){
				for(int dy = -R; dy <= R; dy++){
					if(dx != 0 | dy != 0){
						//neighbor site -> new site = x+dx, y+dy
						int nborSite = ((y+dy+L)%L)*L+((x+dx+L)%L);
						neighborTable[i][nborIndex] = nborSite;
						nborIndex += 1;
					}
				}
			}
		}
	}

}

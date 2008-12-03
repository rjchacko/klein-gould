package rachele.Networks;

import kip.util.Random;
import scikit.jobs.params.Parameters;

public class SmallWorld {
	Random random = new Random();
	public int [] spin;
	public int L;
	int z, R, sumSpins;
	public double T, h, dt, J;
	double time, p;
	int [][] neighborTable;
	
	public SmallWorld(Parameters params) {
		L = Integer.highestOneBit(params.iget("L"));
		params.set("L", L);
		spin = new int [L*L];
		R = params.iget("R");
		z = (2*R+1)*(2*R+1)-1;
		params.set("z",z);
		neighborTable = new int [L*L][z];
		random.setSeed(params.iget("Random seed"));
		setParameters(params);
	}
	
	public void setParameters(Parameters params) {
		T  = params.fget("T");
		J  = params.fget("J", 1);
		h  = params.fget("h", 0.0);
		dt = params.fget("dt");
		p = params.fget("p");
		params.set("magnetization", sumSpins/(double)(L*L));

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
		for(int i = 0; i < L*L*dt; i++){
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
		}
		time += dt;
	}
	
	void findRandomNeigbors(){
		for(int i = 0; i < L*L; i++){
			for(int j = 0; j < z; j++){
				int randomSpin =  (int)(random.nextDouble()*(double)(L*L));
				while (randomSpin == i)
					randomSpin =  (int)(random.nextDouble()*(double)(L*L));
				neighborTable[i][j] = randomSpin;
  			}
		}
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
	
	public void restartClock(){
		time = 0.0;
	}
	
	public void setSquareNeighbors(){	
		System.out.println("Setting Square Neighbors");
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
	
	public void rewire(){
		System.out.println("Rewiring");
		for (int i = 0; i < L*L; i++){
			for (int j = 0; j < z; j++){
				double rand = random.nextDouble();
				if (rand < p){
//					System.out.println(i + " " + j);
					int newNbor = (int)(random.nextDouble()*(double)(L*L));
					while(i == newNbor){
						newNbor = (int)(random.nextDouble()*(double)(L*L));						
					}
					neighborTable[i][j]=newNbor;
				}
			}
		}
	}

}

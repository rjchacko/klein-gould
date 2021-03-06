package rachele.Networks;

//import static java.lang.Math.min;
//import rachele.ising.dim2.IsingLR.DynType;
import scikit.jobs.params.Parameters;
import kip.util.Random;

public class RandomGraph {

	Random random = new Random();
	public int [] spin;
	public int L;
	int z, sumSpins;
	public double T, h, dt, J;
	double time;
	int [][] neighborTable;
	
	public RandomGraph(Parameters params) {
		L = Integer.highestOneBit(params.iget("L"));
		params.set("L", L);
		spin = new int [L*L];
		z = params.iget("z");
		neighborTable = new int [L*L][z];
		random.setSeed(params.iget("Random seed"));
		findRandomNeigbors();
		setParameters(params);
	}
	
	public void setParameters(Parameters params) {
		T  = params.fget("T");
		J  = params.fget("J", 1);
		h  = params.fget("h", 0.0);
		dt = params.fget("dt");
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

}

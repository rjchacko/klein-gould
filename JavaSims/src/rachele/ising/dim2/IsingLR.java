package rachele.ising.dim2;

import static java.lang.Math.exp;
import static java.lang.Math.min;
import kip.ising.RewindableDynamics;
import kip.ising.spinblock.SpinBlocks2D;
import kip.util.Random;
import scikit.jobs.params.Parameters;

public class IsingLR extends RewindableDynamics {
	public SpinBlocks2D spins;


	public enum DynType {METROPOLIS, GLAUBER, KAWA_GLAUBER, KAWA_METROPOLIS};
	public DynType dynamics = DynType.GLAUBER;
	public int L, R;
	public double T, J, h;
	public Random random = new Random();
	public static final double kRpeak = 4.4934092;

	
	public IsingLR(Parameters params) {
		L = Integer.highestOneBit(params.iget("L"));
		params.set("L", L);
		R = min(params.iget("R"), L/2-1);
		params.set("R", R);
		
		spins = new SpinBlocks2D(L, R);
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
		
		dt = params.fget("dt");
		T  = params.fget("T");
		J  = params.fget("J", 1);
		h  = params.fget("h", 0.0);

	}
	
	
	public RewindableDynamics clone() {
		IsingLR ising = (IsingLR)super.clone();
		ising.random = ising.random.clone();
		return ising;
	}
	
	
	public double magnetization() {
		return (double)spins.sumAll() / (L*L);
	}
	

	public void randomizeField(double m) {
		if (m == 1 || m == -1) {
			for (int i = 0; i < L*L; i++)
				spins.set(i%L, i/L, (int)m);
		}
		else {
			for (int i = 0; i < L*L; i++) {
				// p(s = +-1) = (1 +- m) / 2
				int s = (random.nextDouble() < (1+m)/2) ? 1 : -1;
				spins.set(i%L, i/L, s);
			}
		}
	}

	public void restartClock(){
		time = 0.0;
	}
	
	public void setField(double m) {
		double mAcc = 0;
		for (int i = 0; i < L*L; i++) {
			int s = (mAcc > i*m) ? -1 : 1;
			spins.set(i%L, i/L, s);
			mAcc += s;
		}
	}
	
	public double[] getField(int dx) {
		int scale = Integer.numberOfTrailingZeros(dx);
		int blocks[] = spins.blocksAtScale(scale);
		double ret[] = new double[blocks.length];
		double blockSize = (1<<scale)*(1<<scale);
		for (int i = 0; i < blocks.length; i++)
			ret[i] = blocks[i]/blockSize;
		return ret;
	}
	
	
	private boolean shouldFlip(double dE) {
		switch (dynamics) {
			case METROPOLIS:
			case KAWA_METROPOLIS:
				return dE <= 0 || random.nextDouble() < Math.exp(-dE/T);
			case GLAUBER:
			case KAWA_GLAUBER:
				return (dE <= 0 && T == 0) || random.nextDouble() < exp(-dE/T)/(1+exp(-dE/T));
			default:
				assert false;
		}
		return false;
	}
	
	protected void _step() {
		for (int cnt = 0; cnt < L*L*dt; cnt++) {
			int x1 = random.nextInt(L);
			int y1 = random.nextInt(L);
			int s1 = spins.get(x1, y1);
			double dE = 2*s1*(h + J*(spins.sumInRange(x1,y1)-s1)/(4*R*(R+1)));
			switch (dynamics) {
				case METROPOLIS:
				case GLAUBER:
					if (shouldFlip(dE)) {
						spins.flip(x1, y1);
					}
					break;
					
				case KAWA_GLAUBER:
				case KAWA_METROPOLIS:
					int x2 = (x1 + random.nextInt(2*R+1) - R + L)%L;
					int y2 = (y1 + random.nextInt(2*R+1) - R + L)%L;
					int s2 = spins.get(x2, y2);
					if (s2 != s1) {
						dE += 2*s2*(h + J*(spins.sumInRange(x2,y2)-s2)/(4*R*(R+1)));
						if (shouldFlip(dE)) {
							spins.flip(x1, y1);
							spins.flip(x2, y2);
						}
					}
					break;
				default:
					assert false;
			}
			scikit.jobs.Job.yield();
		}
	}
	
	public double dTime(){
		return dt;
	}
	
	public double CountInteractions(){
		double intsPerSpin;
		int intCount = 0;
		for(int x = 0; x < L; x++){
			for (int y = 0; y < L; y ++){
				intCount += (spins.sumInRange(x,y)+(2*R+1)*(2*R+1))/2;
			}
		}
		intCount /= 2;
		intsPerSpin = (double)intCount/(L*L);
		intsPerSpin /=(R*R);
		return intsPerSpin;
	}
	
}

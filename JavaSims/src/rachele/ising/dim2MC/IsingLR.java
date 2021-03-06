package rachele.ising.dim2MC;

import static java.lang.Math.exp;
import static java.lang.Math.min;
import kip.ising.RewindableDynamics;
import kip.ising.spinblock.SpinBlocks2D;
import kip.util.Random;
import scikit.dataset.PointSet;
import scikit.jobs.params.Parameters;

public class IsingLR extends RewindableDynamics {
	public SpinBlocks2D spins;


	public enum DynType {METROPOLIS, GLAUBER, KAWA_GLAUBER, KAWA_METROPOLIS};
	public DynType dynamics = DynType.GLAUBER;
	public int L, R;
	public int jumpRange = 1;
	public int rangeOffset = 0;
	public double T, J, h;
	public Random random = new Random();
	public static final double kRpeak = 4.4934092;

//	public double [] siteEnergy;
	
	public boolean track = false;
	public int noParticles = 0;
	public int [] particleLoc;
	public int [] particleNo;
	
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
		String dyn = params.sget("Dynamics");//, "Ising Glauber");
		if (dyn.equals("Ising Glauber"))
			dynamics = DynType.GLAUBER;
		else if (dyn.equals("Ising Metropolis"))
			dynamics = DynType.METROPOLIS;
		else if (dyn.equals("Kawasaki Glauber")){
			dynamics = DynType.KAWA_GLAUBER;
			jumpRange = params.iget("Jump Range");			
		}else if (dyn.equals("Kawasaki Metropolis")){
			dynamics = DynType.KAWA_METROPOLIS;
			jumpRange = params.iget("Jump Range");			
		}
		
		dt = params.fget("dt");
		T  = params.fget("T");
		J  = params.fget("J", 1);
		h  = params.fget("h", 0.0);

	}
	
	public double get_dt(){
		return dt;
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
				return (dE <= 0 && T == 0.0) || random.nextDouble() < exp(-dE/T)/(1+exp(-dE/T));
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
			double z = (2*R+1)*(2*(R+rangeOffset)+1)-1;

			switch (dynamics) {
				case METROPOLIS:
				case GLAUBER:
					double siteEnergy = -s1*(h + J*(sumWithOffset(x1,y1)-s1)/z);
					double dE = -2*(siteEnergy);
//					if(dE>0 && shouldFlip(dE)) System.out.println("yelp T = " + time());
					if (shouldFlip(dE)) {
						spins.flip(x1, y1);
					}
					break;
					
				case KAWA_GLAUBER:
				case KAWA_METROPOLIS:
					int dx = random.nextInt(2*jumpRange+1) - jumpRange;
					int dy = random.nextInt(2*jumpRange+1) - jumpRange;
					int x2 = (x1 + dx + L)%L;
					int y2 = (y1 + dy + L)%L;
					int s2 = spins.get(x2, y2);
					if (s2 != s1) {
						double siteEnergyBeforeFlip = -s1*(h + J*(sumWithOffset(x1,y1)-spins.get(x1,y1))/z);
						double newSiteEnergyBeforeFlip = -s2*(h + J*(sumWithOffset(x2,y2)-spins.get(x2, y2))/z);
						
						//flip temp
						spins.flip(x1, y1);
						spins.flip(x2, y2);
						
						double siteEnergyAfterFlip = s1*(h + J*(sumWithOffset(x1,y1)-spins.get(x1,y1))/z);
						double newSiteEnergyAfterFlip = s2*(h + J*(sumWithOffset(x2,y2)-spins.get(x2, y2))/z);
						//dE = final energy - initial energy
						dE = (siteEnergyAfterFlip + newSiteEnergyAfterFlip)-(siteEnergyBeforeFlip + newSiteEnergyBeforeFlip);						
//						System.out.println(dE + " should flip " + shouldFlip(dE));
						if (shouldFlip(dE)) {
							if (track) trackParticles(x1,y1,x2,y2);
						}else{
							//flip back and flip site energy
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
	
	public int sumWithOffset(int x, int y) {
		return spins.sumInRange(x-R, x+R, y-R-rangeOffset, y+R+rangeOffset);
	}
	
	public void setRangeOffset(int offset){
		rangeOffset = offset;
	}
	
	private void trackParticles(int x1, int y1, int x2, int y2){
		int site1 = x1+y1*L;
		int site2 = x2+y2*L;
		int oldParticleLoc, newParticleLoc;
		//which site has the particle?
		if(spins.get(x1, y1)==1){
			oldParticleLoc = site1;
			newParticleLoc = site2;
		}else{
			oldParticleLoc = site2;
			newParticleLoc = site1;
		}
		//What is the label of this particle?
		int particleLabel = particleNo[oldParticleLoc];
		if(particleLabel==-1)
			System.out.println("Particle label error");
		//set new particle location
		particleLoc[particleLabel] = newParticleLoc;
		//change labels
		particleNo[oldParticleLoc] = -1;
		particleNo[newParticleLoc] = particleLabel;
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
	
	public PointSet getAveHslice(){
		double slice[] = new double[L];
		for (int x = 0; x < L; x++) {
			int sum=0;
			for(int y = 0; y < L; y++)
				sum += spins.get(x, y);
			slice[x] = sum/(double)L;
		}
		return new PointSet(0, 1.0, slice);
	}
	
	public void initTrackParticles(){
		int pts=0;
		for (int i = 0; i < L*L; i++){
			int x = i%L;
			int y = i/L;
			if (spins.get(x,y)==-1){
				pts += 1;
			}
		}
		noParticles = pts;
		particleNo = new int [L*L];
		particleLoc = new int [noParticles];
		pts = 0;
		for (int i = 0; i < L*L; i++){
			int x = i%L;
			int y = i/L;
			if (spins.get(x,y)==-1){
				particleLoc[pts] = i;
				particleNo[i]=pts;
				pts += 1;
			}else{
				particleNo[i]=-1;
			}
		}
	}
	
}

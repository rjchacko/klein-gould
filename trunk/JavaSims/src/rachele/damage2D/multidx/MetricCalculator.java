package rachele.damage2D.multidx;

public class MetricCalculator {

	public boolean [] alive;
	int pow;
	int L;
	int N;
	double cg_time;
	double lastRecordTime;
	double lastStressTime;
	double [] stressTimeAve;
	double [] dx1StressTimeAve;
	double [] dx1ActTimeAve;
	double [] dx1SizeActTimeAve;
	
	public MetricCalculator(int power, boolean [] liveSites, double tolerance){
		pow = power;
		L = (int) Math.pow(2, pow);
		N = L*L;
		lastRecordTime = 0;
		allocateArrays();
		findAliveArray(liveSites, tolerance);
	}
	


	public void initAveArrays(double d, double [] stress, int [] act, int [] sizeAct){
		for(int i=0; i < N; i++){
			stressTimeAve[i] = stress[i];
			dx1StressTimeAve[i] = stress[i];
			dx1ActTimeAve[i] = (double)act[i];
			dx1SizeActTimeAve[i] = (double)sizeAct[i];
		}
	}
	
	/**
	* Incorporate new stress into average every plate update.
	*/
	public void calcNewStressAveArray(double cg_time, double [] stress){
		double del_t = cg_time -lastStressTime;
		for(int i=0; i < N; i++){
			stressTimeAve[i] = (stressTimeAve[i]*(lastStressTime)+ stress[i]*del_t)/(cg_time);
		}
		lastStressTime = cg_time;
	}
	
	/**
	* Incorporate new stress into average every dt plate updates.
	*/
	public void calcNewDxAveArrays(double cg_time, double [] stress, int [] act, int [] sizeAct){
		double del_t = cg_time -lastRecordTime;
		for(int i=0; i < N; i++){
			dx1StressTimeAve[i] = (dx1StressTimeAve[i]*(lastRecordTime)+ stress[i]*del_t)/(cg_time);
			dx1ActTimeAve[i] = (dx1ActTimeAve[i]*(lastRecordTime) + (double)act[i]*del_t)/cg_time;
			dx1SizeActTimeAve[i] = (dx1SizeActTimeAve[i]*(lastRecordTime) + (double)sizeAct[i]*del_t)/cg_time;
		}
		lastRecordTime = cg_time;
	}
	
	/**
	* Collect the coarse grained (in space and time) metrics.
	* Returns an array with CG stress, activity, and size-activity metric calculations.
	*/
	public double [][] calcCG_metrics(){
		double [][] ret = new double [pow][6];
		for (int i = 0; i < pow; i++){
			ret[i] = dxCalcMets(i);	
		}
		return ret;
	}
	
	public double [] calcStressMets(){
		double [] ret = cgCalcMets(0,stressTimeAve, true);
		return ret;
	}
	
	double [] dxCalcMets(int dxIndex){
		double [] stressRet = cgCalcMets(dxIndex, dx1StressTimeAve, true);
		double [] actRet = cgCalcMets(dxIndex, dx1ActTimeAve, false);
		double [] saRet = cgCalcMets(dxIndex, dx1SizeActTimeAve, false);
		double [] ret = new double [6];
		ret[0] = stressRet[0];
		ret[1] = stressRet[1];
		ret[2] = actRet[0];
		ret[3] = actRet[1];
		ret[4] = saRet[0];
		ret[5] = saRet[1];
		return ret;
	}
	
	double [] cgCalcMets(int index, double [] timeAve, boolean stressArray){
		int dx = (int)Math.pow(2, index);
		double [] g = dxFindArray(dx, timeAve, stressArray);
		int Lp = L/dx;
		int Np = Lp*Lp;
//		System.out.println("index = " + index + " dx = " + dx + " Np = " + Np);
		
		double spaceSum_nd = 0.0;
		int ct = 0;
		for(int i=0; i < Np; i++){
			if(alive[findAliveIndex(i,index)]){
				spaceSum_nd += g[i];
				ct += 1;
			}
		}
		double gbar_nd = spaceSum_nd/(double)ct;
		double metricSum_nd = 0;
		for (int i = 0; i < Np; i++) if(alive[findAliveIndex(i,index)]) metricSum_nd += Math.pow(g[i] - gbar_nd, 2);
		double metric_nd = metricSum_nd/(double)ct;
		
		double spaceSum_wd = 0.0;
		for(int i=0; i < Np; i++) spaceSum_wd += g[i];
		double gbar_wd = spaceSum_wd/Np;
		double metricSum_wd = 0;
		for (int i = 0; i < Np; i++) metricSum_wd += Math.pow(g[i] - gbar_wd, 2);
		double metric_wd = metricSum_wd/Np;
		
		double [] ret = new double [2];
		ret[0] = metric_nd;
		ret[1] = metric_wd;
		return ret;
	}

	double [] dxFindArray(int dx, double [] g, boolean stressArray){
		int Lp=L/dx;
		int Np=Lp*Lp;
		double [] ret = new double [Np];
		for (int i = 0; i<N; i++){
			int cgSite = findCG_site(i, dx);
			ret[cgSite] += g[i];
		}
		// If stress array, calc average.  If activity array, just add activity.
		if(stressArray)
			for (int i = 0; i < Np; i++) ret[i] /= (double)(dx*dx); 
		return ret;
	}
	
	public int findCG_site(int s, int blockSize){
		int x = s % L;
		int y = s / L;
		int xp = x / blockSize;
		int yp = y / blockSize;
		int cgSite = yp *(L/blockSize) + xp;
		return cgSite;
	}
	
	void allocateArrays(){
		stressTimeAve = new double [N];
		dx1StressTimeAve = new double [N];
		dx1ActTimeAve = new double [N];
		dx1SizeActTimeAve = new double [N];

		int al = 0;
		for (int i = 0; i < pow; i++){
			int boxSize = (int) Math.pow(Math.pow(2, i), 2);
			int Np = N/boxSize;
			al += Np;
		}
		alive = new boolean [al];
	}
	
	void findAliveArray(boolean [] liveSites, double tolerance){
		for (int i = 0; i < pow; i++){
			int dx = (int)Math.pow(2, i);
			int boxSize = dx*dx;
			int Np = N/(boxSize);
			int [] aliveCt = new int [Np];
			for(int j = 0; j < N; j++)
				if(liveSites[j]) aliveCt[findCG_site(j,dx)] += 1;
			double [] percentAlive = new double [Np];
			for (int j = 0; j < Np; j++)
				percentAlive[j] = (double)aliveCt[j] / (double)boxSize;
			for (int j = 0; j < Np; j++){
				int index = findAliveIndex(j, i);
				if(percentAlive[j] <= 1.0 - tolerance) alive[index] = false;
				else alive [index] = true;
			}
		}
	}
	
	public int findAliveIndex(int cgSite, int ip){
		int ret = 0;
		for(int i = 0; i < ip; i++){
			int dx = (int)Math.pow(2, i);
			int Np = N/((int) Math.pow(dx, 2));
			ret += Np;
		}
		ret += cgSite;
		return ret;
	}
	
}

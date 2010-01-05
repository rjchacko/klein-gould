package chris.old.ofc;

import chris.util.CopyUtil;
import scikit.jobs.params.Parameters;

public class LowerSr extends Damage2D{

	private double dSr, dSrW;
	private int NfMax, Nfailures[], SiteFailures[], DC;
	private boolean dSrN;
	
	public LowerSr(Parameters params) {
	
		super(params);
		lContructor(params);
		
	}
	
	public void lContructor(Parameters params){
		
		DC = 0;
		dSr = params.fget("d\u03C3_r");
		dSrN = !((dSrW = params.fget("d\u03C3_r width")) == 0);
		return;
	}
	
	public void Initialize(){
		
		super.Initialize();	
		SiteFailures = new int[getN()];
		
	}
	
	protected void InitArrays(){
		
		super.InitArrays();
			
		
		// I guess the failure thresholds are randomized?
		
		double SfNOW = getSf0();
		double SrNOW = getSr0();
		double LoadStress = (1+10E-2)*SrNOW;  // ad hoc!!!!
		
		// reset Sf and stress
		double SrMAX = 0;
		for (int jj = 0 ; jj < getN() ; jj++){
			Sf[jj] = SrNOW + (SfNOW - SrNOW)*rand.nextDouble();
			setStress(jj, LoadStress);	
			if(getSr(jj) > SrMAX) SrMAX = getSr(jj);
			SiteFailures[jj] = 0;
		}
		
		NfMax = (int) (Math.ceil(SrMAX / dSr) + 1);
		if(dSrN){
			if(dSr > dSrW){
				NfMax = (int) (Math.ceil(SrMAX / (dSr-dSrW)) + 1);
			}
			else{
				NfMax = getN();
			}
		}
		else{
			Nfailures = new int[NfMax];
		}
		for (int jj = 0 ; jj < NfMax ; jj++){
			Nfailures[jj] = 0;
		}
		
		Nfailures[0] = getN();
				
	}
	
	protected void ForceFailure(){
	
		if(DC >= getN()){
			seeds = null;
			return;
		}
		
		// pick one site, move its failure threshold, check for failure 
		// do it again if no one fails . . .
		
		// leads to a few notions of time
		
		// set dt by cycle in while loop ... but this F-s up metric calc maybe
		
		
		int jj = rand.nextInt(getN());
		Sf[jj] = getSr(jj) + (Sf[jj] - getSr(jj))*rand.nextDouble();
		
		while(Sf[jj] > getStress(jj)){
			jj = rand.nextInt(getN());
			Sf[jj] = getSr(jj) + (Sf[jj] - getSr(jj))*rand.nextDouble();
		}
		
		seeds = CopyUtil.copyArray(jj,1);
		ManageLives(seeds);	
		adjustNF(seeds);
		
		return;
		
	}
	
	protected double nextSr(int site){
				
		double ret;	
		
		if(Sf[site] > dSr) Sf[site] -= dSr;	// weaken the failure threshold
		
		if(dSrN){
			return ((ret = getSr(site) - dSr - dSrW*rand.nextDouble()) > 0) ? ret : 0;
		}
		else{
			return ((ret = getSr(site) - dSr) > 0) ? ret : 0;
		}
	}
	
	protected void adjustNF(int[] sites){
		
		for(int jj = 0 ; jj < sites.length ; jj++){
			
			if(getSr(sites[jj]) == 0){
				killSite(sites[jj]);
				DC++;
			}

			SiteFailures[sites[jj]]++;
			Nfailures[SiteFailures[sites[jj]]]++;
			Nfailures[SiteFailures[sites[jj]]-1]--;
		}
		
		return;
	}
	
	public int approxFL(){
		int ret = 0;
		for (int jj = 0 ; jj < NfMax ; jj++){
			ret += jj*Nfailures[jj];
		}
		
		return (NfMax - ret);
	}
	

}

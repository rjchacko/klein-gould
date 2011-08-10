package chris.ofcdamage;

import java.util.ArrayList;

import chris.util.MathUtil;

import scikit.jobs.params.Parameters;

public class stepdown2Dfast extends damage2Dfast{
	
	private ArrayList<Integer> liveSites;
	private double h, th;
	private boolean stochastic;

	public stepdown2Dfast(Parameters params){
		
		super(params); // call to damage2Dfast's constructor
		contructor_stepdown2Dfast(params);
		return;
	}
	
	public void contructor_stepdown2Dfast(Parameters params){
		
		h          = 0;
		th         = 1.05;
		stochastic = true;
		for (int jj = 0 ; jj < N ; jj++){
			liveSites.add(jj);
		}
		
		return;
	}
	
	// overrides forceZeroVel from ofc2Dfast
	protected void forceZeroVel(int mct, boolean takedata, boolean eqmode){
		
		if(eqmode){	// *not* running in damage mode
			super.forceZeroVel(mct, takedata, eqmode);
			return;
		}
		else{
			if(stochastic){
				int s0;
				int dt2 = 0;
				double dt = 0;
				double tmp;

				while(true){
					s0 = getRand().nextInt(liveSites.size());
					s0 = liveSites.get(s0);
					tmp    = sf[s0];
					sf[s0] = getSr0() + (sf[s0]-getSr0())*(1+h)*getRand().nextDouble();
					dt     += Math.abs(tmp - sf[s0]);
					dt2    ++;
					if(sf[s0] < th*(getSr0()+getDsr())) // the th is AD HOC
						killSite(s0);
					if(sf[s0] <= stress[s0])
						break;
				}
				tmp = 0;
				if (takedata){
					// calculate metric 
					for(int jj = 0 ; jj < N ; jj++){
						sbar[jj] += MathUtil.bool2bin(!failed[jj])*stress[jj];
						tmp   += MathUtil.bool2bin(!failed[jj])*sbar[jj];
					}
					tmp   = tmp / (N-Ndead);
					Omega = 0;
					for(int jj = 0 ; jj < N ; jj++){
						Omega += MathUtil.bool2bin(!failed[jj])*(sbar[jj] - tmp)*(sbar[jj] - tmp);
					}
					Omega = Omega/((double)(mct)*(double)(mct)*(double)(N)); // multiple dt's !!!!
					saveData(mct, eqmode, dt, dt2);
				}
				fs[newindex++] = s0;
				failSite(s0,mct);	
				return;
			}
			else{
				double dsf, tmpbar;
				int jjmax, ndt;
				index = 0;
				ndt   = 0;
				newindex = index;

				jjmax = 0;
				if(takedata){
					tmpbar = 0;
					Omega  = 0;
					for (int jj = 0 ; jj < N ; jj++){ //use this loop to calculate the metric PART 1
						hstrs.accum(stress[jj]);
						// find next site to fail
						if( ((sf[jj]-stress[jj]) < (sf[jjmax]-stress[jjmax])) && !failed[jj] ) jjmax = jj;
						// calculate metric (PART 1)
						sbar[jj] += MathUtil.bool2bin(!failed[jj])*stress[jj];
						tmpbar   += MathUtil.bool2bin(!failed[jj])*sbar[jj];
						// sum the sites from last time
						ndt += MathUtil.bool2bin(ftt[jj]);
					}
					dsf = sf[jjmax]-stress[jjmax];
					tmpbar = tmpbar / (N-Ndead);
					//tmpbar = tmpbar / N;
					Omega  = 0;
					for (int jj = 0 ; jj < N ; jj++){ //use this loop to calculate the metric PART 2
						// add stress to fail site
						sf[jj] -= MathUtil.bool2bin(!failed[jj])*dsf;
						if(sf[jj] < th*(getSr0()+getDsr())) // the th is AD HOC
							killSite(jj);
						//calculate metric (PART 2)
						Omega += MathUtil.bool2bin(!failed[jj])*(sbar[jj] - tmpbar)*(sbar[jj] - tmpbar);
						// reset ftt
						ftt[jj] = false;
					}
					//calculate metric (PART 3)
					//Omega = Omega/((double)(mct)*(double)(mct)*(double)(N-Ndead));
					Omega = Omega/((double)(mct)*(double)(mct)*(double)(N));


					// save and/or write data
					if(mct%dlength == 0 && mct > 0){
						writeData(mct);
					}
					saveData(mct, eqmode, dsf, ndt);
				}
				else{
					for (int jj = 0 ; jj < N ; jj++){
						// find next site to fail
						if( ((sf[jj]-stress[jj]) < (sf[jjmax]-stress[jjmax])) && !failed[jj] ) jjmax = jj;
					}
					// add stress to fail site
					dsf = sf[jjmax]-stress[jjmax];
					for (int jj = 0 ; jj < N ; jj++){
						sf[jj] -= MathUtil.bool2bin(!failed[jj])*dsf;
						if(sf[jj] < th*(getSr0()+getDsr())) // the th is AD HOC
							killSite(jj);
					}
				}

				fs[newindex++] = jjmax;
				failSite(jjmax,mct);	
				return;
				
			}
		}
	}

	protected void killSite(int site){
		
		super.killSite(site);
		liveSites.remove(liveSites.indexOf(site));
		return;
	}
	
	protected double nextSf(int site){
		double nextsf = getSr0() + (sf[site]-getSr0())*(1+h)*getRand().nextDouble();
		if(nextsf < th*(getSr0()+getDsr())) // the th is AD HOC
			killSite(site);
		return nextsf;
	}
	
}
		

//	m = Math.min(Sf[trialsite]+h*(Sf[trialsite] - Sr[trialsite]),Sf0); 	// h = 0 max is Sf[trialsite] (orig. model in PRE)
//	// h = 1, interval symmetric about Sf[trialsite]
//	// h > 1, healing more probable
//	// h < 1, damage more probable
//trialMove = Sf[trialsite] - (Sr[trialsite] + (m - Sr[trialsite])*rand.nextDouble()); 
	
	
//	int[] temp = new int[getN()];
//	int tempC = 0;
//	
//	if(getTime(-1) < 0){	// deal with equilibration
//		seeds = null;
//		return;		// SfNOW -= dSfNOW;	
//	}
//	
//	SfNOW -= dSfNOW;	
//	
//	for (int jj = 0 ; jj < getN() ; jj++){
//		Sf[jj] = getSr0() + (Sf[jj] - getSr0())*rand.nextDouble();
//		if (getStress(jj) > Sf[jj]){
//			temp[tempC++] = jj;
//		}
//	}
//	
//	seeds = CopyUtil.copyArray(temp,tempC);
//	
//	ManageLives(seeds);
//	setT1T2();	// for purposes of measuring the metric
//	
//	return;
	
	
	////////////////////////////////////////////

//	double dsigma, tmpbar;
//	int jjmax, ndt;
//	index = 0;
//	ndt   = 0;
//	newindex = index;
//
//	jjmax = 0;
//	if(takedata){
//		tmpbar = 0;
//		Omega  = 0;
//		for (int jj = 0 ; jj < N ; jj++){ //use this loop to calculate the metric PART 1
//			hstrs.accum(stress[jj]);
//			// find next site to fail
//			if( ((sf[jj]-stress[jj]) < (sf[jjmax]-stress[jjmax])) && !failed[jj] ) jjmax = jj;
//			// calculate metric (PART 1)
//			sbar[jj] += MathUtil.bool2bin(!failed[jj])*stress[jj];
//			tmpbar   += MathUtil.bool2bin(!failed[jj])*sbar[jj];
//			// sum the sites from last time
//			ndt += MathUtil.bool2bin(ftt[jj]);
//		}
//		dsigma = sf[jjmax]-stress[jjmax];
//		tmpbar = tmpbar / (N-Ndead);
//		//tmpbar = tmpbar / N;
//		Omega  = 0;
//		for (int jj = 0 ; jj < N ; jj++){ //use this loop to calculate the metric PART 2
//			// add stress to fail site
//			stress[jj] += MathUtil.bool2bin(!failed[jj])*dsigma;
//			//calculate metric (PART 2)
//			Omega += MathUtil.bool2bin(!failed[jj])*(sbar[jj] - tmpbar)*(sbar[jj] - tmpbar);
//			// reset ftt
//			ftt[jj] = false;
//		}
//		//calculate metric (PART 3)
//		//Omega = Omega/((double)(mct)*(double)(mct)*(double)(N-Ndead));
//		Omega = Omega/((double)(mct)*(double)(mct)*(double)(N));
//
//		
//		// save and/or write data
//		if(mct%dlength == 0 && mct > 0){
//			writeData(mct);
//		}
//		saveData(mct, eqmode, dsigma, ndt);
//	}
//	else{
//		for (int jj = 0 ; jj < N ; jj++){
//			// find next site to fail
//			if( ((sf[jj]-stress[jj]) < (sf[jjmax]-stress[jjmax])) && !failed[jj] ) jjmax = jj;
//		}
//		// add stress to fail site
//		dsigma = sf[jjmax]-stress[jjmax];
//		for (int jj = 0 ; jj < N ; jj++){
//			stress[jj] += MathUtil.bool2bin(!failed[jj])*dsigma;
//		}
//	}
//	
//	fs[newindex++] = jjmax;
//	failSite(jjmax,mct);	
//	
//	if(catalogue == null)
//		return;
//	catalogue[0][mct] = dsigma;
//	catalogue[1][mct] = jjmax%L;
//	catalogue[2][mct] = (int)(jjmax/L);
//	return;
	


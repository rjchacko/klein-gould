package chris.ofcdamage;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;

import scikit.jobs.params.Parameters;
import chris.util.LatticeNeighbors;
import chris.util.MathUtil;

public class ofc2DfastCG extends ofc2Dfast{

	protected double OmegaCGs,OmegaCGsact, OmegaCGact, sbarCG[], stressCG[], dataCG[][];
	protected int LCG, M, NCG, nbaCG[], act[], actbar[], cgID[], qCG, tau, ssold, sactCG[], sactbarCG[];
	public static int dcatCG = 3;
	
	
	public ofc2DfastCG(Parameters params) {
		
		super(params);
		constructor_ofc2DfastCG(params);
		return;
	}
	
	public void constructor_ofc2DfastCG(Parameters params){
	
		LatticeNeighbors nbsCG;
		
		LCG       = params.iget("L_cg (must be odd)");
		tau       = params.iget("t_cg");
		M         = getL()/LCG;
		NCG       = M*M;
		cgID      = new int[N];
		sbarCG    = new double[NCG];
		stressCG  = new double[NCG];
		sactbarCG = new int[NCG];
		sactCG    = new int[NCG];
		act       = new int[NCG];
		actbar    = new int[NCG];
		nbsCG     = new LatticeNeighbors(getL(),getL(),0,(LCG-1)/2,LatticeNeighbors.Type.PERIODIC,LatticeNeighbors.Shape.Square);
		nbaCG     = nbsCG.get(0);
		qCG       = nbaCG.length;
		OmegaCGs  = 0;
		OmegaCGsact = 0;
		dataCG    = new double[dcatCG][dlength];

		for (int kk = 0 ; kk < M ; kk++){
			for (int jj = 0 ; jj < M ; jj++){
				int site = (LCG-1)/2+jj*LCG+kk*(LCG*getL()); 	// these are the centers of
				cgID[site] = jj+kk*M;
				for (int ii = 0 ; ii < qCG ; ii++){	
					cgID[nbsCG.getJ(site, 0, nbaCG, ii)] = jj+kk*M;
				}
			}
		}

		return;
	}
	
	protected void forceZeroVel(int mct, boolean takedata, boolean eqmode){
		// force failure in the zero velocity limit

		double dsigma, tmpbar, dsmin, tmpbarcgs, tmpbaract, tmpbarsact;
		index    = 0;
		newindex = index;
		dsmin    = 1e5;
		
		if(takedata){
			tmpbar      = 0;
			Omega       = 0;
			tmpbarcgs   = 0;
			tmpbarsact  = 0;
			tmpbaract   = 0;
			OmegaCGs    = 0;
			OmegaCGsact = 0;
			OmegaCGact  = 0;
			
			// add the number of failed sites from last update
			// to size activity metric (the "parent" site is
			// stored in ssold)
			sactCG[cgID[ssold]] += GR;
			ssold = 0;
			for (int jj = 0 ; jj < N ; jj++){ //use this loop to calculate the metric PART 1
				// find next site to fail
				if( (sf[jj]-stress[jj]) < dsmin && !failed[jj]){
					ssold = jj;
					dsmin = sf[jj]-stress[jj];
				}
				// calculate metric (PART 1)
				sbar[jj]           += MathUtil.bool2bin(!failed[jj])*stress[jj];
				tmpbar             += MathUtil.bool2bin(!failed[jj])*sbar[jj];
				stressCG[cgID[jj]] += MathUtil.bool2bin(!failed[jj])*stress[jj];
				if(mct%tau == 0){
					if(jj < NCG){
						sactbarCG[jj] += sactCG[jj];
						tmpbarsact    += sactbarCG[jj];
						actbar[jj]    += act[jj];
						tmpbaract     += actbar[jj];
					}
				}
			}
			dsigma = sf[ssold]-stress[ssold];
			tmpbar = tmpbar / (N-Ndead);
			for (int jj = 0 ; jj < N ; jj++){ //use this loop to calculate the metric PART 2
				// add stress to fail site
				stress[jj] += MathUtil.bool2bin(!failed[jj])*dsigma;
				
				//calculate metric (PART 2)
				Omega += MathUtil.bool2bin(!failed[jj])*(sbar[jj] - tmpbar)*(sbar[jj] - tmpbar);
				if(mct % tau == 0){
					if(jj < NCG){
						if(jj == 0){
							tmpbarsact = (double)(tmpbarsact)/(double)(NCG);
							tmpbaract  = (double)(tmpbaract)/(double)(NCG);
						}
						sbarCG[jj]  += stressCG[jj];
						tmpbarcgs   += sbarCG[jj];
						OmegaCGsact += (tmpbarsact-sactbarCG[jj])*(tmpbarsact-sactbarCG[jj]);
						sactCG[jj]   = 0; //reset counting
						OmegaCGact  += (tmpbaract-actbar[jj])*(tmpbaract-actbar[jj]);
						act[jj]      = 0; //reset counting
					}
					else if(jj < 2*NCG){
						if(jj == NCG) tmpbarcgs = tmpbarcgs/NCG;
						OmegaCGs        += (sbarCG[jj%NCG]-tmpbarcgs)*(sbarCG[jj%NCG]-tmpbarcgs);
						stressCG[jj%NCG] = 0; //reset counting
					}
				}
			}
			//calculate metric (PART 3)
			Omega = Omega/((double)(mct)*(double)(mct)*(double)(N-Ndead));
			if(mct%tau == 0){
				OmegaCGsact  = OmegaCGsact/((double)(mct/tau)*(double)(mct/tau)*(double)(NCG));
				OmegaCGact   = OmegaCGact/((double)(mct/tau)*(double)(mct/tau)*(double)(NCG));
				OmegaCGs     = OmegaCGs/((double)(mct)*(double)(mct)*(double)(NCG));
			}

			// save and/or write data
			if(mct%dlength == 0 && mct > 0){
				writeData(mct);
				writeCGdata(mct);
			}
			act[cgID[ssold]]++;	//activity metric counting
			saveData(mct, eqmode, dsigma);
			saveCGdata(mct);
		}
		else{
			ssold = 0;
			for (int jj = 0 ; jj < N ; jj++){
				// find next site to fail
				if( ((sf[jj]-stress[jj]) < (sf[ssold]-stress[ssold])) && !failed[jj] ) ssold = jj;
			}
			// add stress to fail site
			dsigma = sf[ssold]-stress[ssold];
			for (int jj = 0 ; jj < N ; jj++){
				stress[jj] += MathUtil.bool2bin(!failed[jj])*dsigma;
			}
		}
		
		fs[newindex++] = ssold;
		failSite(ssold,mct);		
		return;
	}
	
	protected void failSite(int index){

		super.failSite(index);
		sactCG[cgID[index]]++;
		return;
	}
	
	protected void failSite(int index, int mct){
		
		failSite(index);
		return;
	}

	protected void saveCGdata(int now){
		
		if(now%tau==0){
			dataCG[0][now%dlength] = 1./OmegaCGs;
			dataCG[1][now%dlength] = 1./OmegaCGsact;
			dataCG[2][now%dlength] = 1./OmegaCGact;
		}
		else{
			dataCG[0][now%dlength] = -1;
			dataCG[1][now%dlength] = -1;
			dataCG[2][now%dlength] = -1;
		}
		return;
	}
	
	public void writeCGdata(int now){
		int ub;
		int offset = (int)((now-1)/dlength);
		offset = offset*dlength;

		ub = now%dlength;
		if(ub==0) ub = dlength;
		
		try{
			File file = new File(getOutdir()+File.separator+getBname()+"_CG.txt");
			PrintWriter pw = new PrintWriter(new FileWriter(file, true), true);
			for (int jj = 0 ; jj < ub ; jj++){
				pw.print(jj+offset);
				pw.print("\t");
				for (int kk = 0 ; kk < dcatCG ; kk++){			
					pw.print(dataCG[kk][jj]);
					pw.print("\t");
				}
				pw.println();
			}			
			pw.close();
		}
		catch (IOException ex){
			ex.printStackTrace();
		}
		return;
	}
	
	public void clearData(){ 
		
		if (N > dlength){
			for (int jj = 0 ; jj < N ; jj++){
				if(jj < dlength){
					for (int kk = 0 ; kk < dcat ; kk++){
						data[kk][jj] = 0;
					}
					for (int kk = 0 ; kk < dcatCG ; kk++){
						dataCG[kk][jj] = 0;
					}		
				}
				sbar[jj] = 0;
				if(jj < NCG){
					sbarCG[jj]    = 0;
					sactbarCG[jj] = 0;
					actbar[jj]    = 0;
					sactCG[jj]    = 0;
					act[jj]       = 0;
				}
			}	
		}
		else{
			for (int jj = 0 ; jj < dlength ; jj++){
				for (int kk = 0 ; kk < dcat ; kk++){
					data[kk][jj] = 0;
				}
				for (int kk = 0 ; kk < dcatCG ; kk++){
					dataCG[kk][jj] = 0;
				}	
				if(jj < N) sbar[jj] = 0;
				if(jj < NCG){
					sbarCG[jj]    = 0;
					sactbarCG[jj] = 0;
					actbar[jj]    = 0;
					sactCG[jj]    = 0;
					act[jj]       = 0;
				}
			}	
		}

		return;
	}
	
	protected void configOF(){
		
		super.configOF();
		try{
			File file = new File(getOutdir()+File.separator+getBname()+"_CG.txt");
			PrintWriter pw = new PrintWriter(new FileWriter(file, true), true);
			pw.print("Time");
			pw.print("\t");
			pw.print("Inverse CG Stress Metric");
			pw.print("\t");
			pw.print("Inverse CG Size Activity Metric");
			pw.print("\t");
			pw.print("Inverse Activity Metric");
			pw.println();
			pw.close();
		}
		catch (IOException ex){
			ex.printStackTrace();
		}
		
		return;
	}
	
}

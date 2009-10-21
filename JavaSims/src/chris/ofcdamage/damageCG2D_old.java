package chris.ofcdamage;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;

import scikit.jobs.params.Parameters;
import chris.util.LatticeNeighbors;
import chris.util.MathUtil;

public class damageCG2D_old extends damage2Dfast{

	protected double OmegaCGs, OmegaCGgr, sbarCG[], stressCG[], grbarCG[], dataCG[][];
	protected int LCG, M, NCG, nbaCG[], grCG[], cgID[], qCG, tau;
	protected LatticeNeighbors nbsCG;
	public static int dcatCG = 2;

	public damageCG2D_old(Parameters params) {

		super(params);
		damageCG2D_old_constructor(params);
		return;
	}

	public void damageCG2D_old_constructor(Parameters params){
	
		LCG       = params.iget("L'");
		tau       = params.iget("\u03C4");
		M         = getL()/LCG;
		NCG       = M*M;
		cgID      = new int[N];
		sbarCG    = new double[NCG];
		stressCG  = new double[NCG];
		grbarCG   = new double[NCG];
		grCG      = new int[NCG];
		nbsCG     = new LatticeNeighbors(getL(),getL(),0,(LCG-1)/2,LatticeNeighbors.Type.PERIODIC,LatticeNeighbors.Shape.Square);
		nbaCG     = nbsCG.get(0);
		qCG       = nbaCG.length;
		OmegaCGs  = 0;
		OmegaCGgr = 0;
		dataCG    = new double[dcatCG][dlength];
		
		for (int kk = 0 ; kk < M ; kk++){
			for (int jj = 0 ; jj < M ; jj++){
				int site = (LCG-1)/2+jj*LCG+kk*(LCG*getL()); 	// these are the centers of
																// of the coarse grained blocks
				for (int ii = 0 ; ii < qCG ; ii++){	
					int nb = getNJ(site, ii);
					if(ii == -1){
						cgID[site] = jj+kk*M;
					}
					else{
						cgID[nb] = jj+kk*M;
					}
				}
			}
		}

		return;
	}
	
	protected void forceFailure(int mct, boolean takedata){
		// force failure in the zero velocity limit
		//

		double dsigma, tmpbar, dsmin, tmpbarcgs;
		int jjmax, tmpbarcgf, ndt;

		index = 0;
		newindex = index;
		dsmin    = 1e5;
		ndt      = 0;
		
		// force failure in the zero velocity limit
		jjmax = 0;
		if(takedata){
			tmpbar    = 0;
			tmpbarcgs = 0;
			tmpbarcgf = 0;
			Omega     = 0;
			OmegaCGs  = 0;
			OmegaCGgr = 0;
			for (int jj = 0 ; jj < N ; jj++){ //use this loop to calculate the metric PART 1
				// find next site to fail (NB, must be alive)
				if( (sf[jj]-stress[jj]) < dsmin && !failed[jj]){
					jjmax = jj;
					dsmin = sf[jj]-stress[jj];
				}
				// calculate metric (PART 1)
				sbar[jj] += MathUtil.bool2bin(!failed[jj])*stress[jj];
				tmpbar   += MathUtil.bool2bin(!failed[jj])*sbar[jj];
				stressCG[cgID[jj]] += MathUtil.bool2bin(!failed[jj])*stress[jj];
				ndt += MathUtil.bool2bin(ftt[jj]);
				if(mct%tau == 0){
					if(jj < NCG){
						grbarCG[jj] += grCG[jj];
						tmpbarcgf   += grbarCG[jj];
					}
				}
			}
			dsigma = sf[jjmax]-stress[jjmax];
			tmpbar = tmpbar / (N-Ndead);
			Omega  = 0;
			for (int jj = 0 ; jj < N ; jj++){ //use this loop to calculate the metric PART 2
				// add stress to fail site
				stress[jj] += MathUtil.bool2bin(!failed[jj])*dsigma;
				//calculate metric (PART 2)
				Omega += MathUtil.bool2bin(!failed[jj])*(sbar[jj] - tmpbar)*(sbar[jj] - tmpbar);
				ftt[jj] = false;
				if(jj < NCG){
					sbarCG[jj]  += stressCG[jj];
					if(mct%tau == 0) OmegaCGgr += (tmpbarcgf-grbarCG[jj])*(tmpbarcgf-grbarCG[jj]);
				}
				else if(jj < 2*NCG){
					tmpbarcgs += sbarCG[jj%NCG];
					if(mct%tau == 0) grCG[jj%NCG] = 0; //reset counting
				}
				else if(jj < 3*NCG){
					OmegaCGs += (sbarCG[jj%NCG]-tmpbarcgs)*(sbarCG[jj%NCG]-tmpbarcgs);
					stressCG[jj%NCG] = 0; // so it starts from zero next time
				}
			}
			//calculate metric (PART 3)
			Omega    = Omega/((double)(mct)*(double)(mct)*(double)(N-Ndead));
			OmegaCGs = OmegaCGs/((double)(mct)*(double)(mct)*(double)(NCG));
			if(mct%tau == 0) OmegaCGgr = OmegaCGgr/((double)(mct/tau)*(double)(mct/tau)*(double)(NCG));

			// save and/or write data
			if(mct%dlength == 0 && mct > 0){
				writeData(mct);
				writeCGdata(mct);
			}
			saveData(mct, false, dsigma, ndt);
			saveCGdata(mct);
		}
		else{
			for (int jj = 0 ; jj < N ; jj++){
				// find next site to fail
				if( (sf[jj]-stress[jj]) < dsmin && !failed[jj]){
					jjmax = jj;
					dsmin = sf[jj]-stress[jj];
				}
			}
			// add stress to fail site
			dsigma = sf[jjmax]-stress[jjmax];
			for (int jj = 0 ; jj < N ; jj++){
				stress[jj] += MathUtil.bool2bin(!failed[jj])*dsigma;
			}
		}
				
		fs[newindex++] = jjmax;
		failSite(jjmax);
		
		return;
	}
	
	protected void forceZeroVel(int mct, boolean takedata){
		
		forceFailure(mct, takedata);
		return;
	}
	
	protected void failSite(int index){

		grCG[cgID[index]]++;
		GR++;
		failed[index]  = true;
		return;
	}
	
	private int getNJ(int site, int nb){

		return nbsCG.getJ(site,0,nbaCG,nb);
	}
	
	protected void saveCGdata(int now){
		
		dataCG[0][now%dlength] = 1./OmegaCGs;
		if(now%tau==0) dataCG[1][now%dlength] = 1./OmegaCGgr; // this HAS to be the last data value
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
				for (int kk = 0 ; kk < (dcatCG-1) ; kk++){			
					pw.print(dataCG[kk][jj]);
					pw.print("\t");
				}
				if(jj%tau==0){
					pw.print(dataCG[dcatCG-1][jj]);
				}
				else{
					pw.print("-");

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
					sbarCG[jj] = 0;
					grbarCG[jj]  = 0;
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
					sbarCG[jj] = 0;
					grbarCG[jj]  = 0;
				}
			}	
		}

		return;
	}
	
	protected void configOF(){
		
		try{
			File file = new File(getOutdir()+File.separator+getBname()+".txt");
			PrintWriter pw = new PrintWriter(new FileWriter(file, true), true);
			pw.print("Time");
			pw.print("\t");
			pw.print("Inverse Metric");
			pw.print("\t");
			pw.print("Shower Size");
			pw.print("\t");
			pw.print("EQ Mode");
			pw.print("\t");
			pw.print("Phi");
			pw.println();
			pw.close();
		}
		catch (IOException ex){
			ex.printStackTrace();
		}
		try{
			File file = new File(getOutdir()+File.separator+getBname()+"_CG.txt");
			PrintWriter pw = new PrintWriter(new FileWriter(file, true), true);
			pw.print("Time");
			pw.print("\t");
			pw.print("Inverse CG Stress Metric");
			pw.print("\t");
			pw.print("Inverse CG Events Metric");
			pw.println();
			pw.close();
		}
		catch (IOException ex){
			ex.printStackTrace();
		}
		
		
		
		return;
	}
}

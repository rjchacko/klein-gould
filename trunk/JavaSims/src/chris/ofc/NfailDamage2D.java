package chris.ofc;


import java.io.File;

import scikit.jobs.Job;
import scikit.jobs.params.Parameters;
import chris.util.DirUtil;
import chris.util.LatticeNeighbors;

public class NfailDamage2D extends SimpleDamage2D{
	
	// Parameters
	 public double Sr0, Sc0, Srwidth, Scwidth, alphawidth, lifewidth, DaboutFS;
	 public int Nlives, rmin, hammersize, Nshowers, SonFSindex;
	 public String lifeshape, residualnoise, criticalnoise, outdir, outfile, PicDir;
	 public double Sr[], Sc[], SsoFar[];
	 public Boolean ShowGrid;
	 public double SonFS[];

	 //String outfileCHK;
	 
	// Constructor
	public NfailDamage2D(Parameters params) {
		
		super(params); // I assume this is a call to SimpleDamage2D's constructor
	
		Sr0           = params.fget("Residual Stress (\u03C3_r)");
		Srwidth       = params.fget("\u03C3_r width");
		Sc0           = params.fget("Critical Stress (\u03C3_c)");
		Scwidth       = params.fget("\u03C3_c width");
		lifewidth     = params.fget("Nlives Width");
		Nlives        = params.iget("Number of Lives");
		rmin          = params.iget("Minimum Interaction Radius (r)");
		hammersize    = params.iget("Hammer Size");		
		lifeshape     = params.sget("Life Style");
		residualnoise = params.sget("\u03C3_r Noise");
		criticalnoise = params.sget("\u03C3_r Noise");
		outdir        = params.sget("Data Directory");

		outfile=outdir+File.separator+"Damage.txt";
		//outfileCHK=outdir+File.separator+"DamageCHK.txt";
		
		PicDir=outdir+"/Pics/";
		DirUtil.MkDir(PicDir);
		
		alphawidth = Nwidth;
		
		alive  = new int[2*N];
		Sr     = new double[N];
		Sc     = new double[N];
		SsoFar = new double[N];
		SonFS  = new double[Nlives*N];
		
		if (Nlives == 1){
			Sr0=0.;
			Srwidth=0.;
		}
		
	}
	
	public void Initialize(String str){
		
		
		
		if (BCs.equals("Bordered")){
			if(shape.equals("Circle")){
				neighbors = new LatticeNeighbors(L,L,rmin,R,LatticeNeighbors.Type.BORDERED,LatticeNeighbors.Shape.Circle);
			}
			else if(shape.equals("Square")){
				neighbors = new LatticeNeighbors(L,L,rmin,R,LatticeNeighbors.Type.BORDERED,LatticeNeighbors.Shape.Square);
			}
			else{
				neighbors = new LatticeNeighbors(L,L,rmin,R,LatticeNeighbors.Type.BORDERED,LatticeNeighbors.Shape.Diamond);
			}
		}
		else{
			if(shape.equals("Circle")){
				neighbors = new LatticeNeighbors(L,L,rmin,R,LatticeNeighbors.Type.PERIODIC,LatticeNeighbors.Shape.Circle);
			}
			else if(shape.equals("Square")){
				neighbors = new LatticeNeighbors(L,L,rmin,R,LatticeNeighbors.Type.PERIODIC,LatticeNeighbors.Shape.Square);
			}
			else{
				neighbors = new LatticeNeighbors(L,L,rmin,R,LatticeNeighbors.Type.PERIODIC,LatticeNeighbors.Shape.Diamond);
			}
		}
				
		if(str.equals("Hammer Blow")){

			LatticeNeighbors Hneighbors;
			
			for (int i = 0 ; i < N ; i++){
				stress[i]  = 0.1*Sc0*rand.nextDouble();
				SsoFar[i]=stress[i];
				alive[i+N] = 1;
				
				if(criticalnoise.equals("On")) {
					Sc[i]=Scwidth*rand.nextGaussian()+Sc0;
				}
				else{
					Sc[i]=Sc0;
				}
				
				
				if (residualnoise.equals("On")) {
					Sr[i]=Srwidth*rand.nextGaussian()+Sr0;
				}
				else{
					Sr[i]=Sr0;
				}
				
				if (lifeshape.equals("Flat")){
					alive[i]=rand.nextInt(Nlives)+1;
				}
				else if (lifeshape.equals("Gaussian")){
					alive[i]=(int)(Math.floor((lifewidth*rand.nextGaussian()+1.*Nlives)+0.5));
				}
				else if(lifeshape.equals("Constant")){
					alive[i]=Nlives;
				}
			}
			
			
			if (L%2 == 1){
				imax = N/2 + 1;
			}
			else{
				imax = N/2 - L/2;
			}
			
			if(hammersize > 1){
				
				if (BCs.equals("Bordered")){
					if(shape.equals("Circle")){
						Hneighbors = new LatticeNeighbors(L,L,0,hammersize,LatticeNeighbors.Type.BORDERED,LatticeNeighbors.Shape.Circle);
					}
					else if(shape.equals("Square")){
						Hneighbors = new LatticeNeighbors(L,L,0,hammersize,LatticeNeighbors.Type.BORDERED,LatticeNeighbors.Shape.Square);
					}
					else{
						Hneighbors = new LatticeNeighbors(L,L,0,hammersize,LatticeNeighbors.Type.BORDERED,LatticeNeighbors.Shape.Diamond);
					}
				}
				else{
					if(shape.equals("Circle")){
						Hneighbors = new LatticeNeighbors(L,L,0,hammersize,LatticeNeighbors.Type.PERIODIC,LatticeNeighbors.Shape.Circle);
					}
					else if(shape.equals("Square")){
						Hneighbors = new LatticeNeighbors(L,L,0,hammersize,LatticeNeighbors.Type.PERIODIC,LatticeNeighbors.Shape.Square);
					}
					else{
						Hneighbors = new LatticeNeighbors(L,L,0,hammersize,LatticeNeighbors.Type.PERIODIC,LatticeNeighbors.Shape.Diamond);
					}
				}
				
				int[] Hnbs = Hneighbors.get(imax);
	
	
				for (int i = 0 ; i < Hnbs.length ; i++){
					stress[Hnbs[i]] += 0.85*Sc[imax];
					SsoFar[Hnbs[i]]=stress[Hnbs[i]];
				}
				
			}
			
		}
		else if(str.equals("Flat")){
			
			for (int i = 0 ; i < N ; i++){
				
				if(criticalnoise.equals("On")) {
					Sc[i]=Scwidth*rand.nextGaussian()+Sc0;
				}
				else{
					Sc[i]=Sc0;
				}
				
				stress[i]  = Sc0*rand.nextDouble();
				if((Sc[i]-stress[i])<(Sc[imax]-stress[imax])) imax=i;
				alive[i+N] = 1;
				
				if (residualnoise.equals("On")) {
					Sr[i]=Srwidth*rand.nextGaussian()+Sr0;
				}
				else{
					Sr[i]=Sr0;
				}
				
				if (lifeshape.equals("Flat")){
					alive[i]=rand.nextInt(Nlives)+1;
				}
				else if (lifeshape.equals("Gaussian")){
					alive[i]=(int)(Math.floor((lifewidth*rand.nextGaussian()+1.*Nlives)+0.5));
				}
				else if(lifeshape.equals("Constant")){
					alive[i]=Nlives;
				}
			}
			
			alive[imax]--;
			alive[imax+N]=0;
			
			// Bring Site with Most Stress to Failure
			
			stressMax = stress[imax];
			for (int i = 0; i<N; i++){
				stress[i]+=Sc[imax]-stressMax;
				SsoFar[i]=stress[i];
			}
			
		}
		else {
			System.out.println("Error! Intialization type " + str + " does not exist!");
		}
	
		
		time=0;
		showernumber=0;
				
		return;
	}
	
	public void Avalanche() {
				
		int Nalive, search;
		SonFSindex=0;
		showernumber=0;
		Nshowers=0;

		// get density of failed sites about imax site
		
		DaboutFS = FSdensity(imax);
		
		time++;

		if(Nbool){

			// redistribute stress from the failed site
						
			int[] nbs = neighbors.get(imax);
	
			Nalive = 0;
			for (int i = 0; i<nbs.length; i++){
				Nalive+=alive[nbs[i]+N];
			}
			if(Nalive>0){
				if(alive[imax]>0){
					release=(stress[imax]-Sr[imax])/Nalive;
				}
				else{
					release=stress[imax]/Nalive;
				}
				for (int i = 0; i<nbs.length; i++){
					stress[nbs[i]]+=(1-alphawidth*rand.nextGaussian()-alpha)*release*alive[nbs[i]+N];
				}
			}
			
			// reset plate conditionally
						
			if (alive[imax]>0) {
				stress[imax]=Sr[imax];
				alive[imax+N]=1;
			}
			else {
				stress[imax]=-2;
			}
			
			// search for the beginning of an avalanche
			
			search=0;
			
			for (int i = 0; i<N ; i++){
				if(stress[i]>Sc[i]){
					SonFS[SonFSindex++]=stress[i];
					dead[search++]=i;
					alive[i+N]=0;
					alive[i]--;
				}
			}
			
			Nshowers+=search;
			
			while (search>0){
				showernumber++;
				// Redistribute avalanche stress
				
				for (int i = 0; i<search; i++){
								
					nbs = neighbors.get(dead[i]);
					Nalive=0;
					for (int j = 0; j<nbs.length; j++){
						Nalive+=alive[nbs[j]+N];
					}			
					if(Nalive>0){									
						if(alive[dead[i]]>0){
							release=(stress[dead[i]]-Sr[dead[i]])/Nalive;
						}
						else{
							release=stress[dead[i]]/Nalive;
						}
						for (int j = 0; j<nbs.length; j++){
							stress[nbs[j]]+=(1-alphawidth*rand.nextGaussian()-alpha)*release*alive[nbs[j]+N];
						}
					}
				}
				
				// reset plate conditionally
				
				for (int i = 0; i<search; i++){
					if (alive[dead[i]]>0) {
						stress[dead[i]]=Sr[dead[i]];
						alive[dead[i]+N]=1;
					}
					else {
						stress[dead[i]]=-2;
					}
				}

				// Look for subsequent avalanche
				
				search=0;
				for (int i = 0; i<N ; i++){
					if(stress[i]>Sc[i]){
						SonFS[SonFSindex++]=stress[i];
						dead[search]=i;
						alive[i+N]=0;
						alive[i]--;
						search++;
					}
				}

				Nshowers+=search;
				showering=true;
				Job.animate();	
			}
			
			
		}
		else{

			// redistribute stress from the failed site
			
			int[] nbs = neighbors.get(imax);
	
			Nalive = 0;
			for (int i = 0; i<nbs.length; i++){
				Nalive+=alive[nbs[i]+N];
			}
			if(Nalive>0){
				if (alive[imax]>0){
					release=(1-alpha)*(stress[imax]-Sr[imax])/Nalive;
				}
				else{
					release=(1-alpha)*stress[imax]/Nalive;
				}
				for (int i = 0; i<nbs.length; i++){
					stress[nbs[i]]+=release*alive[nbs[i]+N];
				}
			}
			
			// reset plate conditionally
			
			if (alive[imax]>0) {
				stress[imax]=Sr[imax];
				alive[imax+N]=1;
			}
			else {
				stress[imax]=-2;
			}
			
			
			// search for the beginning of an avalanche
			
			search=0;
			for (int i = 0; i<N ; i++){
				if(stress[i]>Sc[i]){
					SonFS[SonFSindex++]=stress[i];
					dead[search++]=i;
					alive[i+N]=0;
					alive[i]--;
				}
			}
			
			Nshowers+=search;
			
			while (search>0){
				showernumber++;
				// Redistribute avalanche stress
				
				for (int i = 0; i<search; i++){
					
					nbs = neighbors.get(dead[i]);
					Nalive=0;
					for (int j = 0; j<nbs.length; j++){
						Nalive+=alive[nbs[j]+N];
					}			
					if(Nalive>0){			
						if(alive[dead[i]]>0){
							release=(1-alpha)*(stress[dead[i]]-Sr[dead[i]])/Nalive;
						}
						else{
							release=(1-alpha)*stress[dead[i]]/Nalive;
						}
						for (int j = 0; j<nbs.length; j++){
							stress[nbs[j]]+=release*alive[nbs[j]+N];
						}
					}
				}
				
				// reset plate conditionally
				
				for (int i = 0; i<search; i++){
					if (alive[dead[i]]>0) {
						stress[dead[i]]=Sr[dead[i]];
						alive[dead[i]+N]=1;
					}
					else {
						stress[dead[i]]=-2;
					}
				}

				// Look for subsequent avalanche
				
				search=0;
				for (int i = 0; i<N ; i++){
					if(stress[i]>Sc[i]){
						SonFS[SonFSindex++]=stress[i];
						dead[search]=i;
						alive[i+N]=0;
						alive[i]--;
						search++;
					}
				}
				
				Nshowers+=search;
				showering=true;
				Job.animate();	
			}
			
		}
				
		imax=0;
		for (int i = 0 ; i < N ; i++){
			if((Sc[i]-stress[i])<(Sc[imax]-stress[imax])) imax=i;
		}
		
		if(alive[imax] == 0){
			System.out.println("All Sites Failed!");
			crack=true;
			return;
		}
		
		alive[imax]--;
		alive[imax+N]=0;
		
		stressMax = stress[imax];
		for (int i = 0; i<N; i++){
			stress[i]+=Sc[imax]-stressMax;
		}
		
		showering=false;
		Job.animate();
		
		return;
	}
	
	public int[] LiveSites(){
		
		int[] tmp, rtn;
		int j=0, nonzero=0;
		
		tmp = new int[N];
		
		for( int i = 0 ; i < N ; i++) {
			tmp[i] = alive[i+N];
			if (alive[i+N] != 0) nonzero++;
		}
		
		rtn = new int[nonzero];
		
		for ( int i = 0 ; i < N ; i++){
			if(tmp[i] != 0){
				rtn[j++]=i;
			}
		}
		
		return rtn;
	}
	
	public double GetMax(double[] array){
		
		double max=array[0];
		
		for (int i = 0 ; i < array.length ; i++){
			if(array[i]>max) max=array[i];
		}
		
		
		return max;
	}
	
	public int GetMax(int[] array){
		
		int max=array[0];
		
		for (int i = 0 ; i < array.length ; i++){
			if(array[i]>max) max=array[i];
		}
		
		
		return max;
	}
	
	public double radiusGyration(int index){
		
		long x0,y0,x,y,dx,dy;
		long mass=0;
		double rg=0;
		
		int[] RGnbs = neighbors.get(index);
		
		x0=index%L;
		y0=(int)(index/L);
		if (BCs.equals("Bordered")){
			for (int i=0 ; i<RGnbs.length ; i++){
				x=RGnbs[i]%L;
				y=(int)(RGnbs[i]/L);
				mass+=alive[RGnbs[i]+N];
				dx=(x-x0);
				dy=(y-y0);
				rg+=alive[RGnbs[i]+N]*(dx*dx+dy*dy);
			}
		}
		else{
			for (int i=0 ; i<RGnbs.length ; i++){
				x=RGnbs[i]%L;
				y=(int)(RGnbs[i]/L);
				mass+=alive[RGnbs[i]+N];
				dx=(x-x0);
				dy=(y-y0);
				if (2*Math.abs(dx) > L) dx=L-Math.abs(dx);
				if (2*Math.abs(dy) > L) dy=L-Math.abs(dy);
				rg+=alive[RGnbs[i]+N]*(dx*dx+dy*dy);
			}

		}
		
		if(mass>0){
			rg=rg/mass;
		}
		rg=Math.sqrt(rg);
		
		return rg;
	}
	
	public double EFmetric(){
		
		double ret=0;	
		double Sbar=0;
		
		
		for (int i = 0 ; i < N ; i++){
			SsoFar[i]+=stress[i];
		}
		
		for (int i = 0 ; i < N ; i++){
			Sbar+=SsoFar[i]*alive[i+N];
		}
		
		Sbar=Sbar/N;
		
		for (int i = 0 ; i < N ; i++){
			ret+=(SsoFar[i]*alive[i+N]-Sbar)*(SsoFar[i]*alive[i+N]-Sbar);
		}
		
		ret=ret/(time*time);
		
		return ret;
	}
	
	public double FSdensity(int index){
		
		double ret=0;
		int[] DSnbs = neighbors.get(index);
		
		if (DSnbs.length > 1){
			
			for (int i = 0 ; i < DSnbs.length ; i++){
				ret+=alive[DSnbs[i]+N];
			}
			
			// remove failure of imax 
			ret--;
			
			// divide by all sites save imax
			ret=ret/(DSnbs.length-1);
			
			// ret is the density of life site so use 
			// rho_live + rho_dead = 1
			
			ret=1-ret;

		}
		return ret;
	}
	
	public double GetAve(double set[]){

		double ret=0;
		
		if(set.length > 0){
			
			for (int i = 0 ; i < set.length ; i++){
				ret+=set[i];
			}
			
		} ret=ret/set.length;
		
		return ret;
	}
	
	public double GetAve(double set[], int Lset){

		double ret=0;
		
		if (Lset > 0){
		
			for (int i = 0 ; i < Lset ; i++){
				ret+=set[i];
			}
			
			ret=ret/Lset;
		
		}
		
		return ret;
	}
	
	
}
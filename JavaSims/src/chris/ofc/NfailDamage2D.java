package chris.ofc;


import java.io.File;
import java.io.IOException;
import java.text.DecimalFormat;

import javax.imageio.ImageIO;

import scikit.graphics.dim2.Grid;
import scikit.jobs.Job;
import scikit.jobs.params.Parameters;
import chris.util.DirUtil;
import chris.util.LatticeNeighbors;
import chris.util.PrintUtil;

public class NfailDamage2D extends SimpleDamage2D{
	
	// Parameters
	 public double Sr0, Sc0, Srwidth, Scwidth, alphawidth, lifewidth, DaboutFS;
	 public int Nlives, rmin, hammersize, Nshowers, SonFSindex, NdeadS, Nalive, search;
	 public String lifeshape, residualnoise, criticalnoise, outdir, outfile, PicDir;
	 public double Sr[], Sc[], SsoFar[];
	 public Boolean ShowGrid;
	 public double SonFS[];
	 
	 // Formats
	 
	 public DecimalFormat fmt = new DecimalFormat("0000000");
	 public DecimalFormat fmts = new DecimalFormat("0000");

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
		//hammersize    = params.iget("Hammer Size");		
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
				
		if(str.equals("Flat")){
			
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
		
		SonFSindex=0;
		showernumber=0;
		Nshowers=0;
		
		String PLACEHOLDER = "foobar";
		
		time++;
		
		DaboutFS = FSdensity(imax);
		
		// Distribute stress from initial failure
		
		dead[0]=imax;
		search = 1;
		DistStress(PLACEHOLDER);
		
		// reset plate
		resetPlate();
		
		// search for an avalanche
		findAvalanche();
		
		while (search>0){
			
			showernumber++;
			// redistribute stress of failure(s)
			DistStress(PLACEHOLDER);
			
			// reset plates
			resetPlate();
			
			// look for subsequent avalanche
			findAvalanche();
			
			showering=true;
			Job.animate();
	
		}
		// set up next failure

		
		// Find most stressed site
		imax=0;
		for (int i = 0 ; i < N ; i++){
			if((Sc[i]-stress[i])<(Sc[imax]-stress[imax])) imax=i;
		}
		
		// check if most stressed site has already failed
		if(alive[imax] == 0){
			System.out.println("All Sites Failed!");
			crack=true;
			return;
		}
		
		// kill the site 
		alive[imax]--;
		alive[imax+N]=0;
		
		// bring the most stressed site to failure
		stressMax = stress[imax];
		for (int i = 0; i<N; i++){
			stress[i]+=Sc[imax]-stressMax;
		}
		
		showering=false;
		Job.animate();
				
		return;
	}

	public void DistStress(String HowToDumpStress){

		// HowToDumpStress = evenly
		// angular bias (and then which angular bias)
		//				elipse
		//				uniform + random angle ~ cos\theta
		//				all at one cos\theta
		//				multiple cos\theta # ~ \pi R
		
		
		// if string means to distribute stress evenly . . . 
		
		for (int i = 0; i<search; i++){
			
			int[] nbs = neighbors.get(dead[i]);
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
					
					if (Nbool){
						stress[nbs[j]]+=(1-alphawidth*rand.nextGaussian()-alpha)*release*alive[nbs[j]+N];
					}
					else{
						stress[nbs[j]]+=(1-alpha)*release*alive[nbs[j]+N];
					}
					
				}
			}
		}
		
		return;
	}
	
	public void resetPlate(){

		for (int i = 0; i<search; i++){
			if (alive[dead[i]]>0) {
				stress[dead[i]]=Sr[dead[i]];
				alive[dead[i]+N]=1;
			}
			else {
				stress[dead[i]]=-2;
			}
		}
		
		return;
	}
	
	public void findAvalanche(){
		
		// can really reduce the search to the neighbor sites
		
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
		
		return;
	}

	public int[] LiveSites(){
		
		NdeadS = 0;
		
		int[] tmp, rtn;
		int j=0, nonzero=0;
		
		tmp = new int[N];
		
		for( int i = 0 ; i < N ; i++) {
			tmp[i] = alive[i+N];
			if (alive[i+N] != 0) nonzero++;
			if (alive[i] == 0) NdeadS++;
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
	
	public void TakeData(){
		 double rgyr;
		
		int[] LS = LiveSites(); 
		int LSlength = LS.length;
		if(LSlength>0){
			rgyr=radiusGyration(LS[rand.nextInt(LSlength)]);
		}
		else{
			rgyr=0;
		}
		
		PrintUtil.printlnToFile(outfile,time,Nshowers,NdeadS,rgyr,EFmetric(),GetAve(SonFS,SonFSindex),DaboutFS);
				
		return;
	}
	
	public void TakePicture(Grid grid){
		
		String SaveAs = PicDir + File.separator + grid.getTitle()+fmt.format(time)+"-"+fmts.format(showernumber)+".png";
		
		try {
			ImageIO.write(grid.getImage(), "png", new File(SaveAs));
		} catch (IOException e) {
			System.err.println("Error in Writing File" + SaveAs);
		}
		
		return;
	}

		

	
	
}
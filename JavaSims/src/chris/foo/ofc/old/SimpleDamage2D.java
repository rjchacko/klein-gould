package chris.foo.ofc.old;

import java.util.Random;

import scikit.dataset.Histogram;
import scikit.jobs.Job;
import scikit.jobs.params.Parameters;
import chris.util.LatticeNeighbors;



public class SimpleDamage2D {
	
	// Grid Parameters
	public int L, N;
	public LatticeNeighbors neighbors;
	public Boolean crack;
	
	// Stress Parameters
	public int imax;
	public int alive[], dead[];
	public double alpha, Sc, Sr, release, stressMax;
	public double stress[];
	
	// Interaction Parameters
	public int R, rmin, dnt, noisey;
	public double Nwidth, time;
	public String shape, BCs;
	public boolean Nbool;
	
	// Avalanche Parameters
	public int showernumber;
	public Boolean showering;
	
	// Random Number Generation
	public Random rand;					
	
	// Data
	public Histogram histNS;
	
	// Constructor
	public SimpleDamage2D(Parameters params) {
		
		PseudoConstructorSD(params);
		
	}
	
	public void PseudoConstructorSD(Parameters params){
		
		rand  = new Random(params.iget("Random Seed"));
		L     = params.iget("Lattice Size");
		N     = L*L;
		alpha = params.fget("Dissipation (\u03B1)");
		Sc    = params.fget("Critical Stress (\u03C3_c)");
		R     = params.iget("Interaction Radius (R)");
		Nwidth= params.fget("\u03B1 Width");
		
		stress = new double[N];
		alive  = new int[N];
		dead   = new int[N];
		
		noisey=0;
		Nbool=false;
		if (params.sget("\u03B1 Noise").equals("On")){
			noisey=1;
			Nbool=true;
		}
	
		Nwidth=Nwidth*noisey;

		crack=false;
		
		shape = params.sget("Interaction Shape");
		BCs   = params.sget("Boundary Condtions");
		
		return;
		
	}
	
	
	
	
	// Initialize Stress Lattice
	
	public void Initialize(String str) {
				
		// Initialize All Sites
		
		if (str.equals("Hammer Blow")){
			for (int i = 0; i<N; i++){
				stress[i]=0.1*Sc*rand.nextDouble();	
				alive[i]=1;
			}
			if(N%2==1){
				imax=(int)(N/2)+1;
			}
			else{
				imax=(int)(N/2+L/2);
			}
			stress[imax]=Sc;
			alive[imax]=0;
		}
		else if (str.equals("Flat")){
			imax=0;	
			for (int i = 0; i<N; i++){
				stress[i]=Sc*rand.nextDouble();		// How random??
				if(stress[i]>stress[imax]) imax=i;
				alive[i]=1;
			}
			alive[imax]=0;
		}
		else {
			System.out.println("Error! Intialization type " + str + " does not exist!");
		}
		
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
		
		// Bring Site with Most Stress to Failure
		
		stressMax = stress[imax];
		for (int i = 0; i<N; i++){
			stress[i]+=Sc-stressMax;
		}
		time=0;
		showernumber=0;
	}
	
	public void Avalanche() {
	
		int Nalive;
		int search;
		
		time++;
		
			
		showernumber=0;
			
		// Redistribute Stress
			
		int[] nbs = neighbors.get(imax);
			
		Nalive = 0;
		for (int i = 0; i<nbs.length; i++){
				Nalive+=alive[nbs[i]];
		}
		if(Nalive>0){
			release=stress[imax]/Nalive;
			for (int i = 0; i<nbs.length; i++){
				stress[nbs[i]]+=release*(1-alpha*(Nwidth*rand.nextGaussian()+1.))*alive[nbs[i]];
			}
		}
		stress[imax]=-2;
		
		// Look for avalanche
			
		search=0;
		for (int i = 0; i<N ; i++){
			if(stress[i]>Sc){
				dead[search]=i;
				alive[i]=0;
				search++;
			}
		}
			
		while (search>0) {
			showernumber++;
			showering=true;
			
			// Redistribute avalanche stress
			
			for (int i = 0; i<search; i++){
				
				nbs = neighbors.get(dead[i]);
				Nalive=0;
				for (int j = 0; j<nbs.length; j++){
					Nalive+=alive[nbs[j]];
				}			
				if(Nalive>0){									// what do we do with this stress if "if" fails??
					release=stress[dead[i]]/Nalive;
					for (int j = 0; j<nbs.length; j++){
						stress[nbs[j]]+=release*(1-alpha*(Nwidth*rand.nextGaussian()+1.))*alive[nbs[j]];
					}
				}
				stress[dead[i]]=-2;
			}
			
			// Look for subsequent avalanche
			
			search=0;
			for (int i = 0; i<N ; i++){
				if(stress[i]>Sc){
					dead[search]=i;
					alive[i]=0;
					search++;
				}
			}
		
			histNS.accum((double)(time));
			Job.animate();		
		}

		// Reset stresses to create failure
		
		for (int i = 0; i<N; i++){
			if(stress[i]>=stress[imax]) imax=i;
		}
		
		stressMax = stress[imax];
		
		if(alive[imax]==0){
			System.out.println("All Sites Failed!");
			crack=true;
			return;
		}

		for (int i = 0; i<N; i++){
			stress[i]+=Sc-stressMax;
		}

		alive[imax]=0;
	
		showering=false;		
		
		return;
	}
	
	public int getCenter(int index){
		
		double xc=0;
		double yc=0; 
		int mass=0;
		
		int[] nbs = neighbors.get(index);
		
		for (int i=0 ; i<nbs.length ; i++){
			mass+=alive[nbs[i]];
			xc+=alive[nbs[i]]*nbs[i]%L;
			yc+=alive[nbs[i]]*(int)(nbs[i]/L);
		}
		
		if(mass>0){
			xc=xc/mass;
			yc=yc/mass;
		}

		return (int)(Math.floor(xc+0.5))+((int)(Math.floor(yc+0.5)))*L;
	}
	

}

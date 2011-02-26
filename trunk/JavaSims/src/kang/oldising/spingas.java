package kang.oldising;

import java.awt.Color;

import chris.util.PrintUtil;

import scikit.graphics.ColorPalette;
import scikit.graphics.dim2.Grid;
import scikit.jobs.Control;
import scikit.jobs.Job;
import scikit.jobs.Simulation;
//import scikit.jobs.params.ChoiceValue;
import scikit.jobs.params.DoubleValue;
//import scikit.jobs.params.StringValue;

public class spingas extends Simulation{
	Grid grid1 = new Grid("Spin gas lattice 2d");
	public int L1, L2,M; //parameters for the lattice
	public int i, x, y;  //parameters for the index
	public double T;     //temperature 
	//public double chp;     //chemical potential
	public double K;     //interaction constant after normalization
	public double NK;    //interaction constant before normalization
	
	public int R;   //interaction range
	public int step; //Monte Carlo step
	public int steplimit; //upper limit of MC step
	public int metricstart; //the start of metric calculation
	public int N; //how many points in the metric plot
	public double P; //dilution percentage
	
	public double totalmagnetization;
	public double magnetization;
	public String computer;    //the computer on which i run the program
	
	
	
	public int isingspin[];     //the array of the data
	public int particleposition[]; //tracking the particles
	public int hole[];         // the array for the holes within the range of a particle
	
	
	public double timetotalE[];  
	public double timeaverageE[]; //time average energy for metric
	public double Metric[];  //array for plotting metric
	
	
	public int Nneighber(int a,int i ){// function for the index of nearest neighbor
		int nx,ny; //index for neighbor
		int ni=0;
		nx=(int)i/L2;
		ny=(int)i%L2;
		
		if (a==0) {
			ni=nx*L2+ny-1;
			if  (ny==0) {
				ni=nx*L2+ny+L2-1;
			}
			
		}//(x,y-1) up
		
		if (a==1){
			ni=(nx+1)*L2+ny;
			if  (nx==L1-1) {
				ni=(nx-L1+1)*L2+ny;
			}
			
		}//(x+1,y) right
		
		if (a==2){
			ni=nx*L2+ny+1;
			if  (ny==L2-1) {
				ni=nx*L2+ny-L2+1;
			}
			
		}//(x,y+1) down
		
		if (a==3){
			ni=(nx-1)*L2+ny;
			if  (nx==0) {
				ni=(nx+L1-1)*L2+ny;
			}
		}//(x-1,y) left
		
		return ni;
		
	}
	
	public double interactionE (int j){ //function for interaction energy of nearest neighbors
		double Energy=0;
		int b,k;
		for(b=0; b<4;b++){
			k=Nneighber(b,j);
			if (isingspin[k]!=0)
				Energy=Energy+K*(isingspin[j]+1)*(isingspin[k]+1)/4;
	
		}
		return Energy;	
	}    // this function is not useful
	
	public double longrangeE (int i){
		double Energy=0;
		int totalN=0;
		int nx=i/L2;
		int ny=i%L2;
		int kx, ky;
		
		for (int m=-R; m<=R; m++)
			for (int n=-R; n<=R; n++)
			{
				kx=nx+m;
				ky=ny+n;
				if(nx+m<0)
					kx=nx+m+L1;
				if(nx+m>L1-1)
					kx=nx+m-L1;
				if(ny+n<0)
					ky=ny+n+L2;
				if(ny+n>L2-1)
					ky=ny+n-L2;
				if(isingspin[kx*L2+ky]!=0)
					totalN+=(isingspin[kx*L2+ky]+1)/2;
			}
		Energy=K*(totalN-(isingspin[i]+1)/2)*(isingspin[i]+1)/2;
		return Energy;
	}
	
	
	public static void main (String[] kangliu){
		new Control(new spingas(), "Kang Liu's spin gas model" );
	}
	
	
	public void load(Control kang){
		kang.frame (grid1);
		params.add("lattice's width", 100);
		params.add("lattice's length", 100);
		params.addm("Temperature", new DoubleValue(1, 0, 10).withSlider());
		//params.addm("Chemical potential", new DoubleValue(0, -2, 2).withSlider());
		params.addm("Interaction Constant", new DoubleValue(1, -100,100).withSlider());
		params.add("Interaction range", 10);
		params.add("Monte Carlo step's limit", 1000000);
		params.add("Metric Start at",5000);
		params.add("Metric points", 2000);
		params.add("Diluted Percentage", new DoubleValue(0,0,1).withSlider());
		
		//params.add("Computer", new ChoiceValue("IBM", "Mac"));
		//params.add("Mechanics", new ChoiceValue("Metropolis", "Kawasaki"));
		params.add("MC time");
		params.add("magnetization");
		params.add("Particle number");
		
	}
	
	
    public void animate(){
		
		ColorPalette ising = new ColorPalette ();
		ising.setColor(1, Color.BLACK);   //particles
		ising.setColor(-1, Color.WHITE);    //holes
		ising.setColor(0, Color.RED);
		
		grid1.setColors(ising);
		grid1.registerData(L1, L2, isingspin);
		
		
	}
    
    
    public void clear(){
		grid1.clear();
		
	}
	
    public void run(){
    	int i,j,q;
		R = (int)params.fget("Interaction range");
		steplimit = (int)params.fget("Monte Carlo step's limit");
		metricstart = (int)params.fget("Metric Start at");
		L1 =(int)params.fget("lattice's width");
		L2 =(int)params.fget("lattice's length");
		M = L1 * L2;
		isingspin = new int[M];

		
		

		N = (int)params.fget("Metric points");
		Metric = new double[N];
		
		P = params.fget("Diluted Percentage");
		
		double NMetric=0;
		double totalE=0;
		double averageE=0;// definition for metric variables
		
		int X=0;  //the number of the particles in the spin gas
		
		//randomize the initial state
		for (i=0; i<M; i++){
			double pro = Math.random();
			if (pro < P)
				isingspin[i]=0;               //dilute the lattice first
			if(pro > P)
			{
				isingspin[i]=-1;
				if (Math.random()> 0.5)
				{
				isingspin[i]=1;
				X++;
				}
			}
			
	
		}
		Job.animate();
		params.set("Particle number", X);
		timetotalE = new double[X];
		timeaverageE = new double[X];      //arrays of the energy metric only for the particles
		
		
		// track all the particles
		
		particleposition = new int[X];
		j=0;                     //initialize the index for particleposition array
		for (i=0; i<M; i++){
			if (isingspin[i]==1)
			{
				particleposition[j]=i;         //record the position of the jth particle
				j++;
			}
			
		}
    	
    	//enter the MCS
		
		for (step=0; step< steplimit; step++){
		    T = params.fget("Temperature");
		    NK = params.fget("Interaction Constant");
		    //K=NK/((2*R+1)*(2*R+1)-1);
		    K=NK/(X-1);
		
			for(q=0; q<X; q++)
			{
				int v;  // vth particle
				int k, kx, ky;              //indices for particle position
				int lx, ly;                  //positions for other lattice points within the range
				int h,u;                      //index for hole[] array
				int phole;                   //position of the hole
				int holesnumber;            //the number of the holes in the neighborhood of a particle
				hole = new int [(2*R+1)*(2*R+1)-1] ;  //the array for holes
				v= (int)(Math.random()*X);   //randomly choose a particle
				k= particleposition[v];    //find this particle on the lattice
				kx=(int)k/L2;
				ky=(int)k%L2;               //calculate this particle's position
				
				h=0;
				holesnumber=0;
				for (int m=-R; m<=R; m++)
				{
					for (int n=-R; n<=R; n++)
					{
						
						lx=kx+m;
						ly=ky+n;                              //find the other points
						if(kx+m<0)
							lx=kx+m+L1;
						if(kx+m>L1-1)
							lx=kx+m-L1;
						if(ky+n<0)
							ly=ky+n+L2;
						if(ky+n>L2-1)
							ly=ky+n-L2;                        // check the boundary

						if (isingspin[lx*L2+ly]==-1)
							{
							hole[h]=lx*L2+ly;                  //record the position of the hole
							h++;
							holesnumber++;
							}
						
					}
				}
				if(holesnumber!=0)   // there has to be a hole in the neighborhood if you want to move the particle
				{
					u= (int)(Math.random()*holesnumber);            //randomly choose a hole
					phole= hole[u];                             //find that hole on the lattice
					
					//PrintUtil.printlnToFile("/Users/liukang2002507/Desktop/particle.txt",step, isingspin[k]);
					//if(isingspin[k]==-1)
					//	PrintUtil.printlnToFile("/Users/liukang2002507/Desktop/particle.txt", v,k);
					//PrintUtil.printlnToFile("/Users/liukang2002507/Desktop/hole.txt",step, isingspin[phole]);
					
					double E1=longrangeE(k);
					isingspin[k]=-1;
					isingspin[phole]=1;                //move the particle first and calculate the energy
					double E2=longrangeE(phole);
					
					isingspin[k]=1;
					isingspin[phole]=-1;    //move it back
					
					if (E1>E2)            //if decrease the energy, we should move the particle
					{
						isingspin[k]=-1;
						isingspin[phole]=1;               
						particleposition[v]=phole;         //track the particle
					}
					
					if  (E1<E2)           //if increase the energy
					{
						double probability= Math.random();
						if (probability<=Math.exp((E1-E2)/T))
						{
							isingspin[k]=-1;
							isingspin[phole]=1;          
							particleposition[v]=phole;         //track the particle
						}
						/*if (probability>Math.exp((E1-E2)/T))
						{
							isingspin[k]=1;
							isingspin[phole]=-1;   // in this case, we move the particle back to the original position
						}*/
					}
					
					
				}
				if(step % 10 == 0) params.set("MC time", step);
				Job.animate();


				
			
			}
			
			
			if(step > metricstart){
				
				totalE=0;
				
				for (int y=0; y<X; y++)
				{
				    timetotalE[y]+=longrangeE(particleposition[y]);  // the yth particle's energy at time=step
					
					timeaverageE[y]=timetotalE[y]/(step-metricstart);
				 
				}
				
				for (int x=0; x<X; x++)
				{
						totalE+=timeaverageE[x];
				}
				
				averageE=totalE/X;

				NMetric=0;
				for (int z=0; z< X; z++)
				{
					NMetric+=(timeaverageE[z]-averageE)*(timeaverageE[z]-averageE);
				}
				
			
				if(step<N+metricstart+1)
					{
					Metric[step-metricstart-1]=NMetric/X;
					PrintUtil.printlnToFile("F:/data/spingas/dmetriclarge95.txt",step-metricstart, NMetric/M);
					//PrintUtil.printlnToFile("/Users/liukang2002507/Desktop/metric2.txt",step-metricstart, NMetric/M);
					}
				
				if(step==N+metricstart+1)
				{
					PrintUtil.printlnToFile("F:/data/spingas/dmetric0.txt","");
					PrintUtil.printlnToFile("F:/data/spingas/dmetric0.txt","Lattice length=", L1);
					PrintUtil.printlnToFile("F:/data/spingas/dmetric0.txt","Lattice width=", L2);
					PrintUtil.printlnToFile("F:/data/spingas/dmetric0.txt","K=",NK);
					PrintUtil.printlnToFile("F:/data/spingas/dmetric0.txt","Interaction Range=",R);
					PrintUtil.printlnToFile("F:/data/spingas/dmetric0.txt","Mectric starts at", metricstart);
					PrintUtil.printlnToFile("F:/data/spingas/dmetric0.txt","Final Tempearture=",T);
					PrintUtil.printlnToFile("F:/data/spingas/dmetric0.txt","Diluted Percentage=",P);
					
					/*PrintUtil.printlnToFile("/Users/liukang2002507/Desktop/data/spingas/dmetric1.txt","");
					PrintUtil.printlnToFile("/Users/liukang2002507/Desktop/data/spingas/dmetric1.txt","Lattice length=", L1);
					PrintUtil.printlnToFile("/Users/liukang2002507/Desktop/data/spingas/dmetric1.txt","Lattice width=", L2);
					PrintUtil.printlnToFile("/Users/liukang2002507/Desktop/data/spingas/dmetric1.txt","K=",NK);
					PrintUtil.printlnToFile("/Users/liukang2002507/Desktop/data/spingas/dmetric1.txt","Interaction Range=",R);
					PrintUtil.printlnToFile("/Users/liukang2002507/Desktop/data/spingas/dmetric1.txt","Mectric starts at", metricstart);
					PrintUtil.printlnToFile("/Users/liukang2002507/Desktop/data/spingas/dmetric1.txt","Final Tempearture=",T);
					PrintUtil.printlnToFile("/Users/liukang2002507/Desktop/data/spingas/dmetric1.txt","Diluted Percentage=",P);
					*/
				}
				
				}
			
			
			totalmagnetization=0;	
			for(int s=0; s<M; s++)
			{
				totalmagnetization+=isingspin[s];
			}
			magnetization=totalmagnetization/M;
			params.set("magnetization", magnetization);
			
			
			
			
			}
			
		
		
			
			
			
		}
		
		
		
    	
    }
	

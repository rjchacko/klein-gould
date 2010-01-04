package kang.ising;

import java.awt.Color;
import java.io.File;
import java.io.IOException;
import java.text.DecimalFormat;

import javax.imageio.ImageIO;

import chris.util.PrintUtil;
import chris.util.Random;

//import scikit.graphics.ColorGradient;
import scikit.graphics.ColorPalette;
import scikit.graphics.dim2.Grid;
import scikit.jobs.Control;
import scikit.jobs.Job;
import scikit.jobs.Simulation;
//import scikit.jobs.params.ChoiceValue;
import scikit.jobs.params.DoubleValue;


public class dilutedising extends Simulation{
	
	private DecimalFormat fmt = new DecimalFormat("0000000");

	private Random rand;
	
	Grid grid1 = new Grid("Diluted Ising lattice 2d");
	Grid grid2 = new Grid("Diluted Ising lattice copy");

	public int L1, L2, M; //parameters for the lattice
	public int S1, S2, S; //parameters for the coarsegrain blocks, s is the total number of the spins in a block
	public int B1, B2, B; //parameters for the numbers of the blocks in x and y direction, B is the total number
	public int i, x, y;  //parameters for the index
	public double T;     //temperature 
	public double H;     //field
	public double J;     //interaction constant after normalization
	public double NJ;    //interaction constant before normalization
	public double totalmagnetization;
	public double magnetization;
	public int R;   //interaction range
	public int CR;  //coarse grain range for magnetization
	
	public int deadsites;
	public int step; //Monte Carlo step
	public int steplimit; //upperlimit of MC step
	public int metricstart; //the start of metric calculation
	public int N; //how many points in the metric plot
	public double P;
	
	
	
	public int isingspin[];     //the array of the data
	public int isingcopy[];
	public double coarsegrain[]; //the array of the coarse grain magnetization
	public double CGEnergy[]; // the array of the CGblocks' energy
	
	public double timetotalE[];
	public double timetotalE2[];
	public double timetotalS[];
	
	public double timeaverageE[]; //time average energy for metric
	public double timeaverageE2[];
	public double timeaverageS[];
	public double Metric[];  //array for plotting metric
	public double DoubleMetric[];
	public double SMetric[];  //array for the magnetization metric
	
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
	
	public double interactionE (int j){ //function for interaction energy
		double Energy=0;
		int b,k;
		for(b=0; b<4;b++){
			k=Nneighber(b,j);
			Energy=Energy+J*isingspin[j]*isingspin[k];
		}
		return Energy;
			
	}
	
	public double interactionE2 (int j){ //function for interaction energy
		double Energy=0;
		int b,k;
		for(b=0; b<4;b++){
			k=Nneighber(b,j);
			Energy=Energy+J*isingcopy[j]*isingcopy[k];
		}
		return Energy;
			
	}
	
	
	

	public double coarsegrainE (int i){
		double Energy=0;
		double S=0;
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
				S+=coarsegrain[kx*L2+ky];	
			}
		Energy=J*coarsegrain[i]*S-J*coarsegrain[i]*coarsegrain[i];
		return Energy;
	}
	
	public double longrangeE (int i){
		double Energy=0;
		int S=0;
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
				S+=isingspin[kx*L2+ky];	
			}
		Energy=J*isingspin[i]*S-J;
		return Energy;
	}
	
	public double coarsegrain (int i){
		double phi=0;
		int S=0;
		int nx=i/L2;
		int ny=i%L2;
		int kx, ky;
		
		for (int m=-CR; m<=CR; m++)
			{
			for (int n=-CR; n<=CR; n++)
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
				S+=isingspin[kx*L2+ky];	
			 }
			}
		phi=(double)S/((2*CR+1)*(2*CR+1));
		return phi;
	}          // the function to coarsegrain site[i]
	
	
	public double BlockEnergy (int i){
		int Bx= i/B2;
		int By= i%B2;
		int Sx, Sy;
		int j;
		
		double totalenergy=0;
		double energy=0;
		for(Sx=0;Sx<S1;Sx++)
			for(Sy=0;Sy<S2;Sy++)
			{
				j=(Bx*S1+Sx)*L2+(By*S2+Sy);
				totalenergy+=longrangeE(j);
			}
		energy=totalenergy/S;
		return energy;
		
	}
	
	public static void main (String[] kangliu){
		new Control(new dilutedising(), "Kang Liu's diulted ising model" );
	}
	
	public void load(Control liukang){
		liukang.frame (grid1);
		liukang.frame (grid2);

		params.add("lattice's width", 100);
		params.add("lattice's length", 100);
		params.add("CGBlock's width", 1);
		params.add("CGBlock's length", 1);
		params.addm("Temperature", new DoubleValue(1, 0, 10).withSlider());
		params.addm("Field", new DoubleValue(0, -2, 2).withSlider());
		params.addm("Interaction Constant", 1);
		params.add("Interaction range", 10);
		params.add("Coarsegrain range", 0);

		params.add("Monte Carlo step's limit", 1000000);
		params.add("Metric Start at",5000);
		params.add("Metric points", 2000);
		params.add("Diluted Percentage", new DoubleValue(0,0,1).withSlider());

		params.add("MC time");
		params.add("Metric");
		params.add("DoubleMetric");
		params.add("Magnetization Metric");
		params.add("magnetization");
		
	}
	
	public void animate(){
		
		ColorPalette ising = new ColorPalette ();
		ising.setColor(1, Color.BLACK);
		ising.setColor(-1, Color.WHITE);
		ising.setColor(0, Color.RED);
		
		grid1.setColors(ising);
		grid1.registerData(L1, L2, isingspin);
		grid2.setColors(ising);
		grid2.registerData(L1, L2, isingcopy);
//		String SaveAs = PicDir + File.separator + gridS.getTitle()+fmt.format(time)+".png";

		
		/*String SaveAs = "/Users/cserino/Desktop/p/pic_"+fmt.format(step)+".png";
		try {
			ImageIO.write(grid1.getImage(), "png", new File(SaveAs));
		} catch (IOException e) {
			System.err.println("Error in Writing File" + SaveAs);
		}*/
		
		
		//ColorGradient smooth = new ColorGradient();

		/*for (i = 0; i < M; i++) {
			smooth.getColor(coarsegrain[i], -1., 1.);
		}
		grid2.setColors(smooth);
		grid2.registerData(L1, L1, coarsegrain);
		*/
		
	}
	
	public void clear(){
		grid1.clear();
		grid2.clear();
	
	}
	
	public void run(){
		
		int i,j,k;
		R = (int)params.fget("Interaction range");
		CR = (int)params.fget("Coarsegrain range");
		steplimit = (int)params.fget("Monte Carlo step's limit");
		metricstart = (int)params.fget("Metric Start at");
		
		L1 =(int)params.fget("lattice's width");
		L2 =(int)params.fget("lattice's length");
		M = L1 * L2;
		S1 =(int)params.fget("CGBlock's width");
		S2 =(int)params.fget("CGBlock's width");
		S= S1 * S2;
		B1 = L1/S1;
		B2 = L2/S2;
		B = B1 * B2;

		// rand = new Random(SEED)
//		for(int jj = 0 ; jj < 100 ; jj++)
//			rand.nextDouble();
		//Random newr = rand.clone();
		
		isingspin = new int[M];
	    isingcopy = new int[M];
		CGEnergy = new double[B];   // the energy of the blocks
		
		coarsegrain = new double[M]; // the magnetization of each site after coarsegraining
		
		timetotalE = new double[B];
		timetotalE2 = new double[B];
		timeaverageE = new double[B];
		timeaverageE2 = new double[B];
		timetotalS = new double[B];
		timeaverageS = new double[B];
		

		
		N = (int)params.fget("Metric points");
		Metric = new double[N];
		DoubleMetric = new double[N];
		SMetric = new double[N];
		
		P = params.fget("Diluted Percentage");
		
		double NMetric=0;
		double NMetric2=0;
		double totalE=0;
		double totalE2=0;
		double averageE=0;// definition for metric variables
		
		double NSMetric=0;
		double totalS=0;
		double averageS=0;
		
		deadsites=0;
		
		
		
		//randomize the initial state
		for (i=0; i<M; i++){
			isingspin[i]=-1;
			isingcopy[i]=-1;
	//		if(rand.nextDouble() > 0.5) // in [0,1]
		//		rand.nextInt(25) random integer in [0 24]
			if (Math.random()> 0.5)
				isingspin[i]=1;
			if (Math.random()> 0.5)
				isingcopy[i]=1;
				
		}
		
		for (i=0; i<M; i++){
			if(Math.random()< P)
			{
				isingspin[i]=0;
				isingcopy[i]=0;
				deadsites++;
			}
		}// here, the sequence of dilution and initialization is not important
		
		
		for (step=0; step< steplimit; step++){
			
			for (int f=0; f< M; f++){
			
			    H = params.fget("Field");
			    T = params.fget("Temperature");
			    NJ = params.fget("Interaction Constant");
			
				j=(int) (Math.random()*M); //choose a spin randomly
				

			
				double ZeemanE1=-H*isingspin[j];// initial field energy
				double InterE1=0;
                double InterE2=0;
				
				if (R==0) {
					J=NJ/4;
					InterE1=interactionE(j); // initial interaction energy
					InterE2=-interactionE(j);
				}
				
				if (R!=0) {
					J=NJ/((2*R+1)*(2*R+1)-1);
					InterE1=longrangeE(j);
					InterE2=-longrangeE(j)-2*J;
				}
				
				
				double ZeemanE2=-H*(-isingspin[j]); //field energy after flipping
				
				
				double E1=ZeemanE1+InterE1;
				double E2=ZeemanE2+InterE2;
				
				if (E1>E2){

					isingspin[j] = -isingspin[j];
				}// flip the spin
				
				if (E1<E2){
					if (Math.random()<Math.exp((E1-E2)/T)){

						isingspin[j]=-isingspin[j];
					
					}
				}
				
				//////////////////////////////////////////////////////////////////////////
				
				k=(int) (Math.random()*M);
				
				double ZeemanE21=-H*isingcopy[k];// initial field energy
				double InterE21=0;
                double InterE22=0;
				
				if (R==0) {
					J=NJ/4;
					InterE21=interactionE2(k); // initial interaction energy
					InterE22=-interactionE2(k);
				}
				
				if (R!=0) {
					J=NJ/((2*R+1)*(2*R+1)-1);
					InterE1=longrangeE(k);
					InterE2=-longrangeE(k)-2*J;
				}
				
				
				double ZeemanE22=-H*(-isingspin[k]); //field energy after flipping
				
				
				double E21=ZeemanE21+InterE21;
				double E22=ZeemanE22+InterE22;
				
				if (E21>E22){

					isingcopy[k] = -isingcopy[k];
				}// flip the spin
				
				if (E21<E22){
					if (Math.random()<Math.exp((E21-E22)/T)){

						isingcopy[k]=-isingcopy[k];
					
					}
				}
				
				
				
				
				
			
				if(step % 1 == 0) {
					params.set("MC time", step);
					params.set("Metric", NMetric/(M-deadsites));
					params.set("DoubleMetric", NMetric2/(M-deadsites));
					params.set("Magnetization Metric", NSMetric/(M-deadsites));
				}
				Job.animate();
			}
			
			if(step > metricstart){
				
			if(CR!=0)
			{
				for(int e=0; e<M; e++)
					coarsegrain[e]=coarsegrain(e);
			}
			// first coarsegrain the ising lattice
			
			
			totalE=0;
			totalE2=0;
			totalS=0;
			
			
			
			for (int y=0; y<B; y++)
			{
				if(R==0)
					{
					timetotalE[y]+=interactionE(y);
					timetotalE2[y]+=interactionE2(y);
					timetotalS[y]+=isingspin[y]; 
					}
				if(R!=0)
					{
					if(CR==0)
						{
						if(B==M)
							{
							 timetotalE[y]+=longrangeE(y);
							 timetotalS[y]+=isingspin[y]; 
							}
						
						/*if(B!=M)
							{
							CGEnergy[y]=BlockEnergy(y);
							timetotalE[y]+=CGEnergy[y];
								
							}*/
						
						}
					if(CR!=0)
						timetotalE[y]+=coarsegrainE(y);
					
					}
				
				timeaverageE[y]=timetotalE[y]/(step-metricstart);
				timeaverageE2[y]=timetotalE2[y]/(step-metricstart);
				timeaverageS[y]=timetotalS[y]/(step-metricstart);
			 
			}
			
			
				
		
				for (int x=0; x< B; x++)
			    {

					if(isingspin[x]!=0)
						{
						totalE+=timeaverageE[x];
						totalS+=timeaverageS[x];
						}

			    }
			
			averageE=totalE/(M-deadsites);
			averageS=totalS/(M-deadsites);

			NMetric=0;
			NMetric2=0;
			NSMetric=0;
			
			for (int z=0; z< B; z++)
			{
				if (isingspin[z]!=0)
					{
					
					NMetric+=(timeaverageE[z]-averageE)*(timeaverageE[z]-averageE);
					NMetric2+=(timeaverageE[z]-timeaverageE2[z])*(timeaverageE[z]-timeaverageE2[z]);
					NSMetric+=(timeaverageS[z]-averageS)*(timeaverageS[z]-averageS);
					}
		
				
			}
			
		
			if(step<N+metricstart+1)
				{
				Metric[step-metricstart-1]=NMetric/(M-deadsites);
				DoubleMetric[step-metricstart-1]=NMetric2/(M-deadsites);
				SMetric[step-metricstart-1]=NSMetric/(M-deadsites);


				PrintUtil.printlnToFile("F:/data/dilutedising/e-metric8.txt",step-metricstart, Metric[step-metricstart-1]);
				PrintUtil.printlnToFile("F:/data/dilutedising/d-metric8.txt",step-metricstart, DoubleMetric[step-metricstart-1]);
				PrintUtil.printlnToFile("F:/data/dilutedising/s-metric8.txt",step-metricstart, SMetric[step-metricstart-1]);

				//PrintUtil.printlnToFile("/Users/cserino/Desktop/metric2.txt",step-metricstart, NMetric/B);
				}
			
			
			if(step==metricstart+500)
			{
				PrintUtil.printlnToFile("F:/data/dilutedising/e-metric1.txt","");
				PrintUtil.printlnToFile("F:/data/dilutedising/e-metric1.txt","Lattice length=", L1);
				PrintUtil.printlnToFile("F:/data/dilutedising/e-metric1.txt","Lattice width=", L2);
				PrintUtil.printlnToFile("F:/data/dilutedising/e-metric1.txt","Block length=", S1);
				PrintUtil.printlnToFile("F:/data/dilutedising/e-metric1.txt","Block width=", S2);
				PrintUtil.printlnToFile("F:/data/dilutedising/e-metric1.txt","J=",NJ);
				PrintUtil.printlnToFile("F:/data/dilutedising/e-metric1.txt","Interaction Range=",R);
				PrintUtil.printlnToFile("F:/data/dilutedising/e-metric1.txt","Mectric starts at", metricstart);
				PrintUtil.printlnToFile("F:/data/dilutedising/e-metric1.txt","Final Tempearture=",T);
				PrintUtil.printlnToFile("F:/data/dilutedising/e-metric1.txt","Final Field=",H);
				PrintUtil.printlnToFile("F:/data/dilutedising/e-metric1.txt","Diluted Percentage=",P);
				
				PrintUtil.printlnToFile("F:/data/dilutedising/s-metric1.txt","");
				PrintUtil.printlnToFile("F:/data/dilutedising/s-metric1.txt","Lattice length=", L1);
				PrintUtil.printlnToFile("F:/data/dilutedising/s-metric1.txt","Lattice width=", L2);
				PrintUtil.printlnToFile("F:/data/dilutedising/s-metric1.txt","Block length=", S1);
				PrintUtil.printlnToFile("F:/data/dilutedising/s-metric1.txt","Block width=", S2);
				PrintUtil.printlnToFile("F:/data/dilutedising/s-metric1.txt","J=",NJ);
				PrintUtil.printlnToFile("F:/data/dilutedising/s-metric1.txt","Interaction Range=",R);
				PrintUtil.printlnToFile("F:/data/dilutedising/s-metric1.txt","Mectric starts at", metricstart);
				PrintUtil.printlnToFile("F:/data/dilutedising/s-metric1.txt","Final Tempearture=",T);
				PrintUtil.printlnToFile("F:/data/dilutedising/s-metric1.txt","Final Field=",H);
				PrintUtil.printlnToFile("F:/data/dilutedising/s-metric1.txt","Diluted Percentage=",P);
				
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
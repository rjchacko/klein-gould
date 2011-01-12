package kang.ising;

import java.awt.Color;
import java.io.File;
import java.io.IOException;
import java.text.DecimalFormat;

import javax.imageio.ImageIO;

import chris.util.PrintUtil;
import chris.util.Random;

import scikit.graphics.ColorGradient;
import scikit.graphics.ColorPalette;
import scikit.graphics.dim2.Grid;
import scikit.jobs.Control;
import scikit.jobs.Job;
import scikit.jobs.Simulation;
import scikit.jobs.params.DoubleValue;

import static java.lang.Math.PI;

import kang.ising.StructureFactorOpt;

public class antiferroising extends Simulation{
	
	private DecimalFormat fmt = new DecimalFormat("0000");
	
	Grid grid1 = new Grid("Antiferromagnetic Ising lattice 2d for nucleation");
	Grid grid2 = new Grid("Antiferromagnetic Ising initial copy");
	Grid grid3 = new Grid("critical droplet");
	Grid grid4 = new Grid("intervention");
	Grid grid5 = new Grid("Dilution heatmap");
	Grid grid6 = new Grid("Structure Factor");
	

	

	public int L1, L2, M; //parameters for the lattice
    public int i, x, y;  //parameters for the index
	public double T;     //temperature
	public double QuenchT;  //the temperature after the quench
	public double QuenchH;  
	public double H;     //field
	public double field;
	public double J;     //interaction constant after normalization
	public double NJ;    //interaction constant before normalization
	public double percent;  //dilution percent


	public double magnetization;
	public double Ms;  //saturated magnetization
	public int R;   //interaction range R=0 is NN interaction
	
	//bias variables
	public int biasrange;
	public double biaspercent;
	

	//public int numberofcopies;
	public int lifetime;   // the time starts at the moment when quench of the fields happens
	public int quenchstart; // the time of the quench

	public int deadsites;  //the number of diluted sites in the system
	
	public int prestep;  // step for preparation of intial configuration
	public int step; //Monte Carlo step
	public int rationalMCS; // use this we can run half of the MCS
	public int steplimit; //upperlimit of MC step
	public int meanlifetime;
	public int totallifetime;
	
	public int growthpercentage;
	public int decaynumber;
	public int grownumber;
	public int interventionsteplimit;
	public int interventioncopies;
	public double interventionstart;    // this double format is used for pinpoint the critical droplet
	public double interventionP;
	
	

	public int isingspin[];     //the array of the data
	public int initialcopy[];   //the array of the initial copy of the system
	public int isingcopy[];     //the array for reset
	public int interventioncopy[];
	public double dilutionmap[];  

	public Random spinfliprand;
    public Random spinrand;
	public Random dilutionrand;
	public Random newpositionrand;
	public Random newspinfliprand;
	public Random biasrand;
	

	public int spinflipseed;
	public int spinseed;
	public int dilutionseed;
	public int biasseed;
	public int type;  // initialization type(0-random spin    1-all spins up    2-all spins down)
	
	public StructureFactorOpt SFO;
	public double sf[];
	public double fftspin[];
	public double circularSF[];   // S(q)
	public double squareSF[];
	public int Lp;    //elements per side for Structure factor
	public double L;
	public double dq;
	
	
	// functions for the structure factor
	
	
	public double circularaverage(double SF[], int Lp, double q, double dq)   //returns the circular average of SFactor between q and q+2*dq*q
	{
		double totalSF=0;
		int SFcount=0;
		int ix,iy;
		double qi;
		
		for(int i=0; i<Lp*Lp; i++)
		{
			ix=(i%Lp)-Lp/2;
			iy=(i/Lp)-Lp/2;
			qi=ix*ix+iy*iy;
			if((qi>=q*q)&(qi<=((q+dq)*(q+dq))))
			{
				totalSF+=SF[i];
				SFcount++;
			}
		}
		
		return totalSF/SFcount;
	}
	
	public double squareaverage(double SF[], int Lp, double q) //returns the square average of structure factor between q to q+dq
	{
		double totalSF=0;
		int i1,i2,i3,i4;
		int qi;
		qi=(int)q;
		i1=(qi+Lp/2)+(Lp/2)*Lp;
		i2=(-qi+Lp/2)+(Lp/2)*Lp;
		i3=(Lp/2)+(qi+Lp/2)*Lp;
		i4=(Lp/2)+(-qi+Lp/2)*Lp;
		totalSF=SF[i1]+SF[i2]+SF[i3]+SF[i4];
	
		return totalSF/4;
	}
	
	
	public void circularSF(double SF[], int Lp, double dq)
	{
		int qN=0;
		qN= (int) (Lp/(2*dq));
		circularSF= new double [qN];
		for(int qq=0; qq< qN; qq++)
		{
			circularSF[qq]=circularaverage(SF, Lp, qq*dq, dq);
		}
	}
	
	public void squareSF(double SF[], int Lp)
	{
		int qN=0;
		qN = Lp/2;
		squareSF= new double[qN];
		for(int qq=0;qq<qN; qq++)
		{
			squareSF[qq]=squareaverage(SF, Lp , qq);
		}
	}
	
	public void printcircularSF(double CSF[], int Lp, double dq)
	{
		int qN=0;
		qN= (int) (Lp/(2*dq));
		for(int qq=0; qq< qN; qq++)
		{
			if(CSF[qq]>=0.0001)
				PrintUtil.printlnToFile("/Users/liukang2002507/Desktop/circularSF.txt", qq*dq*(2*PI*R)/L, CSF[qq]);
		}
		
	}
	
	public void printsquareSF(double SSF[], int Lp)
	{
		int qN=0;
		qN=Lp/2;
		for(int qq=0; qq<qN; qq++)
		{
			if(SSF[qq]>=0.0001)
				PrintUtil.printlnToFile("/Users/liukang2002507/Desktop/squareSF.txt", qq*(2*PI*R)/L, SSF[qq]);
				
		}
	}
	
	public void printXSF(double SF[], int Lp)
	{
		for(int qx=0; qx<Lp/2; qx++)
			{
			int ix=(qx+Lp/2);
			PrintUtil.printlnToFile("/Users/liukang2002507/Desktop/XSF.txt", qx, qx*(2*PI*R)/L, SF[ix]);
			}
		
	}
	
	public void printYSF(double SF[], int Lp)
	{
		for(int qy=0; qy<Lp/2; qy++)
			{
			int iy=(qy+Lp/2)*Lp;
			PrintUtil.printlnToFile("/Users/liukang2002507/Desktop/YSF.txt", qy, qy*(2*PI*R)/L, SF[iy]);
			}
		
	}
	
	public void findpeak(double SF[], int Lp)
	{
		double max=0;
		int maxposition=0;
		for(int i=0; i<Lp*Lp; i++)
		{
			if(SF[i]>max)
			{
				max=SF[i];
				maxposition=i;
			}
		}
		int x=maxposition%Lp-Lp/2;
		int y=maxposition/Lp-Lp/2;
		PrintUtil.printlnToFile("/Users/liukang2002507/Desktop/MAX.txt", x, y, SF[maxposition]);
		PrintUtil.printlnToFile("/Users/liukang2002507/Desktop/MAX.txt", x*(2*PI*R)/L, y*(2*PI*R)/L, SF[maxposition]);
		
	}
	
	public void printallSF(double SF[], int Lp)
	{
	    for(int i=0; i<Lp*Lp; i++)
	    {
	    	int x=i%Lp-Lp/2;
	    	int y=i/Lp-Lp/2;
	    	PrintUtil.printlnToFile("/Users/liukang2002507/Desktop/all.txt", x, y, SF[i]);
	    	PrintUtil.printlnToFile("/Users/liukang2002507/Desktop/all.txt", x*(2*PI*R)/L, y*(2*PI*R)/L, SF[i]);
	    }
	}
	
	
	
	
 	public int Nneighber(int a,int i )
 	{// function for the index of nearest neighbor
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

 	public double interactionEchange (int range, double constJ, int spin[], int j)
 	{ //function for interaction energy
		double Energy=0;
		double Energychange=0;
		if (range==0)
		{
			J= constJ/4;                     // normalize the interaction constant
			int b,k;
		    for(b=0; b<4;b++)
		    {
			k=Nneighber(b,j);
			Energy=Energy+J*spin[j]*spin[k];
		    }
		    Energychange=-2*Energy;
		    
		}
		if (range!=0)
		{
			int S=0;
			int nx=j/L2;
			int ny=j%L2;
			int kx, ky;
			J= constJ/((2*range+1)*(2*range+1)-1);
			for (int m=-range; m<=range; m++)
				for (int n=-range; n<=range; n++)
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
					S+=spin[kx*L2+ky];	
				}
			Energy=J*spin[j]*S-J;
			Energychange=-2*Energy;
			
		}
		
		return Energychange;	
    }

 	public void spinflip(int range, double constJ, int spin[], int j, Random flip, double temperature, double field)
	{
		double ZeemanE=2*field*spin[j]; //the change in Zeeman's energy if we flip the spin[j]
		double Echange=ZeemanE+interactionEchange(range, constJ, spin,j);
		
		if(Echange<0)
		{
			spin[j]=-spin[j];
		}
		
		else
		{
			if(flip.nextDouble()<=Math.exp(-Echange/temperature))
					spin[j]=-spin[j];
		}
		
	}

 	public double magnetization(int spin[])   // function to calculate the magnetization for spin[]
	{
		double totalmagnetization=0;	
		for(int s=0; s<M; s++)
		{
			totalmagnetization+=spin[s];
		}
		double m=totalmagnetization/M;
		return m;
	}
	
	public void MCS(int spin[], Random flip, double ratio)
	{
	    int j=0;
		H = params.fget("Field");
	    T = params.fget("Temperature");
	    NJ = params.fget("Interaction Constant");
	    rationalMCS= (int) (ratio*M);
	    for (int f=0; f< rationalMCS; f++)
	    {
		   j=(int) (flip.nextDouble()*M);
		   spinflip(R, NJ, spin, j, flip, T, H);
		  
	    }
	    magnetization=magnetization(spin);
	    params.set("magnetization", magnetization);
	    
	}

	public void temperaturequench(double finaltemp)
	{
		params.set("Temperature", finaltemp);
	}
	public void fieldset(double finalfield)
	{
		params.set("Field", finalfield);
	}
	
	
	public void fieldquench()
	{
	    double finalfield=-H;
	    params.set("Field", finalfield);
	    
	}

	public void animate()
	{
		ColorPalette ising = new ColorPalette ();
		ising.setColor(1, Color.BLACK);
		ising.setColor(-1, Color.WHITE);
		ising.setColor(0, Color.RED);
		ising.setColor(2, Color.BLUE);
		ising.setColor(-2, Color.GREEN);
		
		ColorGradient heatmap = new ColorGradient();
		
		grid1.setColors(ising);
		grid1.registerData(L1, L2, isingspin);
		grid2.setColors(ising);
		grid2.registerData(L1, L2, initialcopy);
		grid3.setColors(ising);
		grid3.registerData(L1, L2, interventioncopy);
		grid4.setColors(ising);
		grid4.registerData(L1, L2, isingcopy);
		grid5.setColors(heatmap);
		grid5.registerData(L1, L2, dilutionmap);
		grid6.setColors(heatmap);
		grid6.registerData(Lp, Lp, sf);

	}

	public void clear()
	{
		grid1.clear();
		grid2.clear();
		grid3.clear();
		grid4.clear();
		grid5.clear();
		grid6.clear();

	}
	
	public void initialize(int spin[], int type, double p, int spinseed, int biasR, int dilutionseed, double biasp, int biasseed)
	{
		spinrand= new Random(spinseed);
		dilutionrand= new Random(dilutionseed);
		biasrand= new Random(biasseed);
		deadsites=0;
		
		int biasrangelabel[];
		biasrangelabel= new int [M];    // the array to label the location of bias diluted site(1 for bias diluted site)
		
		for (int b=0;b<M;b++)
			biasrangelabel[b]=0;
			
		
		int cx, cy; //x,y index for the center
		cx=(int)L1/2;
		cy=(int)L2/2;
		if(biasp!=p)
			//square(biasrangelabel, biasR, cx, cy);
		    rectangle(biasrangelabel, 64, biasR, cx,cy);
		
		for(int t=0; t<M; t++)
			spin[t]=1;
		
		
		for(int j=0; j<M; j++)
		{
			if (biasrangelabel[j]==1)
				if(biasrand.nextDouble()<biasp)
					{
					spin[j]=0;
					deadsites++;
					}
			if (biasrangelabel[j]==0)
				if(dilutionrand.nextDouble()<p)
					{
					spin[j]=0;
					deadsites++;
					}
		}
		
		if(type==0)
			for (i=0; i<M; i++)
		    {
			   if(spin[i]!=0)
			   {
				   spin[i]=-1;
			   if (spinrand.nextDouble()> 0.5)
				   spin[i]=1;
			   }
				
		    }
		if(type==1)
			for (i=0; i<M; i++)
		    {
			   if(spin[i]!=0)
				   spin[i]=1;
		    }
		if(type==2)
			for (i=0; i<M; i++)
		    {
			   if(spin[i]!=0)
				   spin[i]=-1;
		    }
		
		for (i=0; i<M; i++){

			initialcopy[i]=spin[i];  // here, make the copy of the system
			
		}    
		
		
	}

	public int X(int bx)
	{
		int realx=bx;
		if (bx>L1)
			realx=bx-L1;
		if (bx<0)
			realx=bx+L1;
		return realx;
	}
	
	public int Y(int by)
	{
		int realy=by;
		if (by>L2)
			realy=by-L2;
		if (by<0)
			realy=by+L2;
		return realy;
	}

	public void rectangle(int label[], int a, int b, int cx, int cy)
	
	{
		int bx, by;
		int x,y;
		for(int t=0; t<M; t++)
			label[t]=0;
		for(bx=cx-a; bx<cx+a; bx++)
			for(by=cy-b; by<cy+b; by++)
			{
				x=X(bx);
				y=Y(by);
				label[x*L2+y]=1;
			}
		
	}
	
	//public void diagonal(int label[],)
	
	public void square(int label[], int range, int cx, int cy)  //the function to draw a square
	{
		int bx, by;
		int x, y;
		for(int t=0; t<M; t++)
			label[t]=0;
		for(bx=cx-range; bx<cx+range; bx++)
			for(by=cy-range; by<cy+range; by++)
			{
				x=X(bx);
				y=Y(by);
				label[x*L2+y]=1;
			}
		
	}
	
	
	public int distance (int a, int b)     // the code to calculate the distance between two points on the lattice
	{
		int dis=0;
		int ax, ay, bx, by;
		int dx2, dy2;
		ax= a/L2;
		ay= a%L2;
		bx= b/L2;
		by= b%L2;
		
		dx2=(ax-bx)*(ax-bx);
		dy2=(ay-by)*(ay-by);
		if((ax-bx+L1)*(ax-bx+L1)<(ax-bx)*(ax-bx))
			dx2=(ax-bx+L1)*(ax-bx+L1);
		if((ax-bx-L1)*(ax-bx-L1)<(ax-bx)*(ax-bx))
			dx2=(ax-bx-L1)*(ax-bx-L1);
		if((ay-by+L2)*(ay-by+L2)<(ay-by)*(ay-by))
			dy2=(ay-by+L2)*(ay-by+L2);
		if((ay-by-L2)*(ay-by-L2)<(ay-by)*(ay-by))
			dy2=(ay-by-L2)*(ay-by-L2);

		dis=dx2+dy2;
		return dis;
	}

	public void movie(Grid grid, int number, int copynumber)   //function to capture the grid
	{
		
			String SaveAs = "/Users/liukang2002507/Desktop/simulation/antiferro/pic_"+fmt.format(copynumber)+"_"+fmt.format(number)+".png";
		try {
			ImageIO.write(grid.getImage(), "png", new File(SaveAs));
		} catch (IOException e) {
			System.err.println("Error in Writing File" + SaveAs);
		}
		
	}
	
	public void intervention(int copies, int steplimit, double percentage)
	{
		//growthpercentage=0;
		decaynumber=0;
		grownumber=0;
		for(int k=0; k<copies; k++)
		{
			//newpositionrand=new Random(k);
			newspinfliprand=new Random(k+spinflipseed);

			for(int l=0; l<M; l++)
				isingcopy[l]=interventioncopy[l];
			for(int s=0; s<steplimit; s++)
			{
				MCS(isingcopy, newspinfliprand, 1);
				Job.animate();
				
			}
			movie(grid4, steplimit, k);
			if(magnetization>(percentage*Ms))
				decaynumber++;
			else
				grownumber++;
			
			params.set("u#copies", k);
			params.set("grow",grownumber);
			params.set("decay",decaynumber);
			
			
		}
		
		
	}

	public double dilutionratio(int spin[], int r, int i, int total)
	{
		double ratio;
		double dilutedsite;
		dilutedsite=0;
		double totalinrange;
		totalinrange=0;
		int j;
		int disij;

		for(j=0; j<total;j++)
		{
            disij= distance(i,j);

			if(disij<=r*r)
			{
				totalinrange++;
				if(spin[j]==0)
					dilutedsite++;
			}
		}
	
		ratio=dilutedsite/totalinrange;
		return ratio;
	}

	public void dilutionheatmap(int spin[], int r)
	{
		double map[];
		map= new double[M];
		for(int s=0; s<M; s++)
		{
			map[s]=dilutionratio(spin, r, s, M);
			dilutionmap[s]=map[s];
		}
		
	}
	
	public static void main (String[] antiferroising){
		new Control(new antiferroising(), "Kang Liu's antiferromagnetic ising model" );
	}
	
	public void load(Control antiferroising){
		antiferroising.frame (grid1);
		antiferroising.frame (grid2);
		antiferroising.frame (grid3);
		antiferroising.frame (grid4);
		antiferroising.frame (grid5);
		antiferroising.frame (grid6);
		
		
		params.add("lattice's width", 128);
		params.add("lattice's length", 128);
		params.add("Lp", 128);
		params.add("dq", new DoubleValue(0.1,0,50).withSlider());
		params.add("Diluted Percentage", new DoubleValue(0.1,0,1).withSlider());
		params.add("Bias percent", new DoubleValue(0.1, 0, 1).withSlider());
		
		params.add("Quench starts at", 100);
		
		params.addm("Interaction Constant", 1);
		params.add("Interaction range", 46);
		params.add("Bias range", 10);
		params.add("Monte Carlo step's limit", 1000000);
		
		params.addm("Temperature", new DoubleValue(9, 0, 10).withSlider());
		params.addm("Quench temperature", new DoubleValue(0.06, 0, 10).withSlider());
		
		params.addm("Field", new DoubleValue(0.75, -2, 2).withSlider());
		params.addm("Quench field", new DoubleValue(0.75, -2, 2).withSlider());

		
		
		params.add("intervention start", new DoubleValue(2880, 0, 99999999).withSlider());
		params.add("intervention steplimit", 50);
		params.add("intervention copies", 50);
		params.add("intervention percentage", new DoubleValue(0.85,0,1).withSlider());
		
		///below is the parameters that would be displayed on the panel
		
		params.add("MC time");
		params.add("Metastable state time");
		params.add("magnetization");
		params.add("Saturated magnetization");
		params.add("lifetime");

        params.add("u#copies");
        params.add("grow");
        params.add("decay");

		
	}
	
    public void run(){
		
    	//int estart;   //start of the cluster evolution
    	
		
		R = (int)params.fget("Interaction range");
		biasrange= (int)params.fget("Bias range");
	    steplimit = (int)params.fget("Monte Carlo step's limit");
		quenchstart=(int)params.fget("Quench starts at");
		QuenchT=params.fget("Quench temperature");
		QuenchH=params.fget("Quench field");

		L1 =(int)params.fget("lattice's width");
		L2 =(int)params.fget("lattice's length");
		Lp =(int)params.fget("Lp");
		dq = params.fget("dq");
		M = L1 * L2 ;
		L=L1;

		
		isingspin = new int[M];
		initialcopy = new int[M];
		isingcopy = new int[M];
		interventioncopy = new int[M];
		dilutionmap = new double [M];
		sf= new double [Lp*Lp];
		fftspin= new double[M];
		
		SFO=new StructureFactorOpt(Lp,L);

		spinseed=1;
		dilutionseed=1;
		biasseed=1;
	
		percent = params.fget("Diluted Percentage");
		biaspercent = params.fget("Bias percent");
		type=0;
		field = params.fget("Field");
		
		interventionstart=params.fget("intervention start");
		interventioncopies=(int)params.fget("intervention copies");
		interventionsteplimit=(int)params.fget("intervention steplimit");
		interventionP=params.fget("intervention percentage");
		
		initialize(isingspin, type, percent, spinseed, biasrange,dilutionseed, biaspercent, biasseed);

		Job.animate();
		
		spinflipseed = (int)((field+1)*QuenchT*10000);
		spinfliprand= new Random (spinflipseed);
		
		for (prestep=0; prestep < 4900000; prestep++)
		{
		//temperaturequench(9);	
		MCS(isingspin, spinfliprand ,1);
		params.set("MC time", prestep-50);
		Job.animate();
		}
		

		temperaturequench(QuenchT);
		Ms=magnetization;
		fieldset(QuenchH);
		
		
		for(step=0; step<30; step++)
		{
			for(int rs=0; rs<100; rs++)
			{
			MCS(isingspin, spinfliprand, 0.01);
			/*for(int ft=0; ft<M; ft++)
			{
				fftspin[ft]=isingspin[ft];
			}
			SFO.takeFT(fftspin);
			for(int e=0; e<Lp*Lp; e++)
			{
				sf[e]=SFO.sFactor[e];
			}
			*/
			params.set("MC time", step);
			Job.animate();
			if((step==-1)&(rs==111250))
			{
				circularSF(SFO.sFactor,Lp,dq);
				printcircularSF(circularSF,Lp,dq);
				squareSF(SFO.sFactor, Lp);
				printsquareSF(squareSF, Lp);
				//printXSF(SFO.sFactor, Lp);
				//printYSF(SFO.sFactor, Lp);
				//findpeak(SFO.sFactor, Lp);
				//printallSF(SFO.sFactor, Lp);
			}
			
			movie(grid1, rs+1 ,step+1);	
			}
		}
		
		for(int astep=0; astep<steplimit; astep++)
		{
			MCS(isingspin, spinfliprand, 1);
			params.set("MC time", astep+20);
			Job.animate();
			//movie(grid1, 9999,astep);	
			
		}
		
	
	
    }// the end of run()
	
	
	
}

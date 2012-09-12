package kang.AEM.BasicStructure;


import chris.util.Random;

public class AEMStructure{
	
	//structures and parameters
	public double wealth[];
	public double initialcopy[];   //the array of the initial copy of the system
	public double skill[];  //the quantity of the agent's to make money in the trade
	
    public double totalwealth;  //the total wealth of the whole system
	public double meanwealth; //the average wealth of the system
    
	public int L1, L2, M; //parameters for the lattice                                                                                                                                                                                                                                                                                                                                                                                                        

	public double percent;   //the trading percent of the less wealth holder's total wealth
	public double tax;  //the universal tax rate during the trade
	
	public double Ngrowth;   // the amount of the wealth growth for every step
	public double growth;  //the amount of the wealth growth after every trading event
    public double alpha;  // this represent the dissipation for the tax
   
	public int R;   //interaction range R=0 is NN interaction
	
    public double order; // a parameter to measure the inequality

	
	//the function for this class IsingStructure
	
	public AEMStructure(int L1, int L2, int R, double percent, double tax, double alpha, double Ngrowth)     //generating function
	{
		this.L1=L1;
		this.L2=L2;
		this.M=L1*L2;
		this.R=R;
		this.percent=percent;
		this.tax=tax;
		this.alpha=alpha;
		this.Ngrowth=Ngrowth;
		this.growth=Ngrowth/this.M;
		
	
		this.wealth=new double [M];
		this.initialcopy=new double [M];
		this.skill=new double [M];
		
	
		
	}
	
	public AEMStructure clone()
	{
		AEMStructure copy= new AEMStructure(L1, L2, R, percent, tax, alpha, Ngrowth);
		for(int t=0;t<M; t++)
		{
			copy.wealth[t]=wealth[t];
			copy.initialcopy[t]=initialcopy[t];
			copy.skill[t]=skill[t];
			
		}
		copy.totalwealth=totalwealth;
		copy.meanwealth=meanwealth;

		return copy;
	}
	
	public void Rinitialization(int Sseed, double min, double max)//wealth configuation randomly initialization
	{
		
		
		Random rand= new Random(Sseed);
	    totalwealth=0;

		for (int i=0; i<M; i++)
		    {
			   wealth[i]=min+rand.nextDouble()*(max-min);
			   initialcopy[i]=wealth[i];  // here, make the copy of the system
			   totalwealth+=wealth[i];
		    }
	
		
	}
	
	public void Uinitialization(double wi)     //wealth configuration uniformly initialization
	{
		for (int i=0; i<M; i++)
	    {
		   wealth[i]=wi;
		   initialcopy[i]=wealth[i];  // here, make the copy of the system
	    }
		totalwealth=wi*M;
	}
	
	public int X(int bx)
	{
		int realx=bx;
		if (bx>=L1)
			realx=bx-L1;
		if (bx<0)
			realx=bx+L1;
		return realx;
	}
	
	public int Y(int by)
	{
		int realy=by;
		if (by>=L2)
			realy=by-L2;
		if (by<0)
			realy=by+L2;
		return realy;
	}

	public int Nneighber(int a,int i )// function for the index of nearest neighbor
 	{
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
	
	public int findInRange(int j, int R, Random rand)  //randomly find a trading target in the range of (2R+1)*(2R+1) square or (R+1)^2-R^2 diamond
	{
		
		int kx, ky;
		int target;
		
	    if(R==L2)
	    {
	    	target=(j+rand.nextInt(M))%M;             //for fully connected model
	    }
	    else if(R==0)
	    {
	    	target=Nneighber(rand.nextInt(4), j);       //for nearest neighbor model
	    }
	    
	    else                       // for long range square model
	    {
	    	int nx=j/L2;
			int ny=j%L2;
			

	        int rx=rand.nextInt(2*R);
	        int ry=rand.nextInt(2*R);
			
	        while((rx==R)&&(ry==R))
	        {
	        	 rx=rand.nextInt(2*R);
	             ry=rand.nextInt(2*R);
	        }
	        
	        kx=X(nx+rx-R);
			ky=Y(ny+ry-R);
			target=kx*L2+ky;
			
	    }
		
	    return target;

	}
	
	public void TSfast(Random flip, double percent, double tax, double alpha, double Ngrowth)// fast trading step  (growth after each step instead of each transaction)
	{
        int j=0;
        int target=0;
        int richer=0;    // not used, just a label for the richer
        int poorer=0;
        double tradeamount=0;
        double aftertax=0;
        double totaltax=0;
	    double totalorder=0;
	   
	    for (int f=0; f< M; f++)
	    {
		   j=flip.nextInt(M);
		   target=findInRange(j, R, flip);    //choose the trading target, then decide who is richer
		   if(wealth[j]>=wealth[target])
		   {
			   richer=j;
			   poorer=target;
		   }
		   else
		   {
			   richer=target;
			   poorer=j;
		   }
		   //after deciding who is the poorer, we can decide the trading amount=percent*wealth[poorer]
		   tradeamount=percent*wealth[poorer];
		   aftertax=(1-tax)*tradeamount;
		   totaltax+=tax*tradeamount;
		   
		   
		   //now decide who win this trade
		   if(flip.nextBoolean())    //assume this means j lose
		   {
			   wealth[j]-=aftertax;
			   wealth[target]+=aftertax;
		   }
		   else     // this means j wins
		   {
			   wealth[j]+=aftertax;
			   wealth[target]-=aftertax;
		   }
		   
	    }
	    
		for(int g=0; g<M; g++)   //now everyone's wealth will grow with the same amount =growth
		   {
			   wealth[g]+=Ngrowth;
			   
		   }
		
		   totalwealth+=Ngrowth*M;
		   meanwealth+=Ngrowth;
	    
	    
	    for(int ii=0; ii<M; ii++)
	    {
	    	wealth[ii]+=totaltax*(1-alpha)/M;       //after all the trading within one step, everybody gets a benefit from the tax after a dissipation alpha
	    	
	    	totalorder+=Math.log(wealth[ii]/meanwealth);
	    }
	    
	    order=totalorder/M;
	    
	}
	
	
	
	public void TS(Random flip, double percent, double tax, double alpha, double growth)// trading step
	{
        int j=0;
        int target=0;
        int richer=0;    // not used, just a label for the richer
        int poorer=0;
        double tradeamount=0;
        double aftertax=0;
        double totaltax=0;
	    double totalorder=0;
	   
	    for (int f=0; f< M; f++)
	    {
		   j=flip.nextInt(M);
		   target=findInRange(j, R, flip);    //choose the trading target, then decide who is richer
		   if(wealth[j]>=wealth[target])
		   {
			   richer=j;
			   poorer=target;
		   }
		   else
		   {
			   richer=target;
			   poorer=j;
		   }
		   //after deciding who is the poorer, we can decide the trading amount=percent*wealth[poorer]
		   tradeamount=percent*wealth[poorer];
		   aftertax=(1-tax)*tradeamount;
		   totaltax+=tax*tradeamount;
		   
		   
		   //now decide who win this trade
		   if(flip.nextBoolean())    //assume this means j lose
		   {
			   wealth[j]-=aftertax;
			   wealth[target]+=aftertax;
		   }
		   else     // this means j wins
		   {
			   wealth[j]+=aftertax;
			   wealth[target]-=aftertax;
		   }
		   
		   
		   for(int g=0; g<M; g++)   //now everyone's wealth will grow with the same amount =growth
		   {
			   wealth[g]+=growth;
			   
		   }
		   totalwealth+=growth*M;
		   meanwealth+=growth;
	    }
	    
	    
	    
	    for(int ii=0; ii<M; ii++)
	    {
	    	wealth[ii]+=totaltax*(1-alpha)/M;       //after all the trading within one step, everybody gets a benefit from the tax after a dissipation alpha
	    	
	    	totalorder+=Math.log(wealth[ii]/meanwealth);
	    }
	    
	    order=totalorder/M;
	    
	}
	

	
}

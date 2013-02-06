package kang.AEM.BasicStructure;


import chris.util.Random;

public class AEMStructure{
	
	//structures and parameters
	public double wealth[];
	public double initialcopy[];   //the array of the initial copy of the system
	public double skill[];  //the quantity of the agent's to make money in the trade
	public double growthmatrix[];   // the distribution of the capability to grow (growthmatrix[i]~wealth[i]^gamma where gamma is the parameter for the growth mode. gamma=1 is linear with wealth and gamma=0 is uniform growth)
	public double normalization; //sum over the growthmatrix[] as a normalization factor
	
	
    public double totalwealth;  //the total wealth of the whole system
	public double meanwealth; //the average wealth of the system
    
	public int L1, L2, M; //parameters for the lattice                                                                                                                                                                                                                                                                                                                                                                                                        

	public double percent;   //the trading percent of the less wealth holder's total wealth
	public double tax;  //the universal tax rate during the trade
	
	public double Ngrowth;   // the amount of the wealth growth for every step
	public double growth;  //the amount of the wealth growth after every trading event
    public double alpha;  // this represent the dissipation for the tax
   
	public int R;   //interaction range R=0 is NN interaction
	
    public double order; // a parameter to measure the inequality (entropy!=order   entropy=order/N+ln(N))
 
    public double sumtrading; // the total trading amount in a trading step
    public double sumflow;   //the total amount of wealth flow towards the rich people, negative represent the flow is towards the poor
   
	
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
		meanwealth=wi;
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
	
	public void SublinearTS(Random flip, double percent, double tax, double alpha, double Ngrowth, double gamma)  //the rich gets richer in a sublinear way using the growthmatrix[]
	{
        int j=0;
        int target=0;
        int richer=0;    // not used, just a label for the richer
        int poorer=0;
        double tradeamount=0;
        double aftertax=0;
        double totaltax=0;
	    double totalorder=0;
	    growthmatrix=new double[M];
	    
	    sumtrading=0;
	    sumflow=0;
	    int direction=1;  //the parameter to determine the direction of the flow
	   
	    double totalgrowth= Ngrowth*M;
	    
	    
	    for (int f=0; f< M; f++)
	    {
		   j=flip.nextInt(M);
		   target=findInRange(j, R, flip);    //choose the trading target, then decide who is richer
		   if(wealth[j]>=wealth[target])
		   {
			   richer=j;
			   poorer=target;
			   direction=1;
		   }
		   else
		   {
			   richer=target;
			   poorer=j;
			   direction=-1;  //because j is the poor, so if j wins, the flow is negative
		   }
		   //after deciding who is the poorer, we can decide the trading amount=percent*wealth[poorer]
		   tradeamount=percent*wealth[poorer];
		   aftertax=(1-tax)*tradeamount;
		   totaltax+=tax*tradeamount;	 
		   
		   sumtrading+=aftertax;
		   
		   //now decide who win this trade
		   if(flip.nextBoolean())    //assume this means j lose
		   {
			   wealth[j]-=aftertax;
			   wealth[target]+=aftertax;
			   sumflow-=direction*aftertax;
		   }
		   else     // this means j wins
		   {
			   wealth[j]+=aftertax;
			   wealth[target]-=aftertax;
			   sumflow+=direction*aftertax;
		   }
		   
	    }
	    
	    totalwealth=0;
	    
	    
	    //now generate the growthmatrix
	    normalization=0;
	    for(int gm=0; gm<M; gm++)
	    {
	    	growthmatrix[gm]=Math.pow(wealth[gm], gamma);
	    	normalization+=growthmatrix[gm];
	    	
	    }
	    
	    
	    
		for(int g=0; g<M; g++)   //now everyone's wealth will grow with the amount proportional to growthmatrix[]
		   {
			   wealth[g]=wealth[g]+growthmatrix[g]*totalgrowth/normalization;
			   totalwealth+=wealth[g];
			   
		   }
		
	
		   meanwealth=totalwealth/M;
	    
	    
	    for(int ii=0; ii<M; ii++)
	    {
	    	wealth[ii]+=totaltax*(1-alpha)/M;       //after all the trading within one step, everybody gets a benefit from the tax after a dissipation alpha
	    	
	    	if(wealth[ii]!=0)
	    		totalorder-=wealth[ii]/meanwealth*Math.log(wealth[ii]/meanwealth);    //might be problematic if there is tax because of the change in the total wealth
	    }
	    
	    order=totalorder;
	}
	
	public void UnfairTSfast(Random flip, double percent, double tax, double alpha, double pmu)  //the rich gets richer based on the growth mode(wi=wi*(1+mu))
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
	    
	    totalwealth=0;
		for(int g=0; g<M; g++)   //the rich gets richer based on the growth mode(wi=wi*(1+pmu))
		   {
			   wealth[g]=wealth[g]*(1+pmu);
			   totalwealth+=wealth[g];
		   }
		
	
		   meanwealth=totalwealth/M;
	    
	    
	    for(int ii=0; ii<M; ii++)
	    {
	    	wealth[ii]+=totaltax*(1-alpha)/M;       //after all the trading within one step, everybody gets a benefit from the tax after a dissipation alpha
	    	
	    	if(wealth[ii]!=0)
	    		totalorder-=wealth[ii]/meanwealth*Math.log(wealth[ii]/meanwealth);    //might be problematic if there is tax because of the change in the total wealth
	    }
	    
	    order=totalorder;
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
	    	
	    	if(wealth[ii]!=0)
	    		totalorder-=wealth[ii]/meanwealth*Math.log(wealth[ii]/meanwealth);
	    }
	    
	    order=totalorder;
	    
	}
	
	public void BiasTSfast(double biasp, Random flip, double percent, double tax, double alpha, double Ngrowth)// bias toward richer fast trading step  (growth after each step instead of each transaction)
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
		   if(flip.nextDouble()<=biasp)    //assume this means poor lose
		   {
			   wealth[poorer]-=aftertax;
			   wealth[richer]+=aftertax;
		   }
		   else     // this means poor wins
		   {
			   wealth[poorer]+=aftertax;
			   wealth[richer]-=aftertax;
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
	    	if(wealth[ii]!=0)
	    		totalorder-=wealth[ii]/meanwealth*Math.log(wealth[ii]/meanwealth);
	    }
	    
	    order=totalorder;
	    
	}
	
	public void TSfastRandom(Random flip, double min, double max, double tax, double alpha, double Ngrowth)// fast trading step  (growth after each step instead of each transaction)
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
		   tradeamount=flip.nextDouble()*(max-min)+min;
		   aftertax=(1-tax)*tradeamount;
		   totaltax+=tax*tradeamount;
		   
		   if(wealth[poorer]>=tradeamount)
		   {
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
	    	if(wealth[ii]!=0)
	    		totalorder-=wealth[ii]/meanwealth*Math.log(wealth[ii]/meanwealth);
	    }
	    
	    order=totalorder;
	    
	}
	
	public void TSfastFix(Random flip, double fixamount, double tax, double alpha, double Ngrowth)// fast trading step  (growth after each step instead of each transaction)
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
		   tradeamount=fixamount;
		   aftertax=(1-tax)*tradeamount;
		   totaltax+=tax*tradeamount;
		   
		   if(wealth[poorer]>=tradeamount)
		   {
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
	    	
	    	if(wealth[ii]!=0)
	    		totalorder-=wealth[ii]/meanwealth*Math.log(wealth[ii]/meanwealth);
	    }
	    
	    order=totalorder;
	    
	}
	
	public void BiasTSfastFix(double biasp, Random flip, double fixamount, double tax, double alpha, double Ngrowth)// fast trading step  (growth after each step instead of each transaction)
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
		   tradeamount=fixamount;
		   aftertax=(1-tax)*tradeamount;
		   totaltax+=tax*tradeamount;
		   
		   if(wealth[poorer]>=tradeamount)
		   {
			 //now decide who win this trade
			   if(flip.nextDouble()<=biasp)    //assume this means poor lose
			   {
				   wealth[poorer]-=aftertax;
				   wealth[richer]+=aftertax;
			   }
			   else     // this means poor wins
			   {
				   wealth[poorer]+=aftertax;
				   wealth[richer]-=aftertax;
			   }
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
	    	
	    	if(wealth[ii]!=0)
	    		totalorder-=wealth[ii]/meanwealth*Math.log(wealth[ii]/meanwealth);
	    }
	    
	    order=totalorder;
	    
	}
	
	public void TSfastFixExp(Random flip, double fixamountpercent, double tax, double alpha, double Ngrowth)// fast trading step  (growth after each step instead of each transaction)
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
		   tradeamount=fixamountpercent*meanwealth;
		   aftertax=(1-tax)*tradeamount;
		   totaltax+=tax*tradeamount;
		   
		   if(wealth[poorer]>=tradeamount)
		   {
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
	    	
	    	if(wealth[ii]!=0)
	    		totalorder-=wealth[ii]/meanwealth*Math.log(wealth[ii]/meanwealth);
	    }
	    
	    order=totalorder;
	    
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
	    	
	    	if(wealth[ii]!=0)
	    		totalorder-=wealth[ii]/meanwealth*Math.log(wealth[ii]/meanwealth);
	    }
	    
	    order=totalorder;
	    
	}
	
	
	public double[][] findrich(double data[], int number)   // the function to find the richest agents in the array
	{
		double[][] richest=new double[2][number];  //richest[0][] are the index of the rich agents and richest[1][] are their wealth
		double threshold=data[0];
		
		int doorman=0;     //the index for the poorest in the rich ones
		for(int i=0; i<number; i++)
		{
			richest[0][i]=i;
			richest[1][i]=data[i];
			if(data[i]<threshold)    
				{
				threshold=data[i];
				doorman=i;                  //find the doorman, the poorest in the first n agents
				}
		}
		
		for(int j=number; j<data.length; j++)
		{
			
			if(data[j]>threshold)
			{
				richest[1][doorman]=data[j];    //replace the doorman with the new member
				richest[0][doorman]=j;
				//now find the new doorman
				threshold=richest[1][0];
				doorman=0;
				for(int k=0; k<number; k++)
				{
					if(richest[1][k]<threshold)
						{
						threshold=richest[1][k];
						doorman=k;
						}
				}
			}
			
		}
		
		return richest;
	}
	
	public double[][] findpoor(double data[], int number)   // the function to find the poorest agents in the array
	{
		double[][] poorest=new double[2][number];
		double threshold=data[0];
		
		int doorman=0;     //the index for the poorest in the rich ones
		for(int i=0; i<number; i++)
		{
			poorest[0][i]=i;
			poorest[1][i]=data[i];
			if(data[i]>threshold)
				{
				threshold=data[i];
				doorman=i;                  //find the doorman
				}
		}
		for(int j=number; j<data.length; j++)
		{
			
			if(data[j]<threshold)
			{
				poorest[1][doorman]=data[j];    //replace the doorman with the new member
				poorest[0][doorman]=j;
				
				//now find the new doorman
				threshold=poorest[1][0];
				doorman=0;
				for(int k=0; k<number; k++)
				{
					if(poorest[1][k]>threshold)
						{
						threshold=poorest[1][k];
						doorman=k;
						}
				}
			}
			
		}
		
		return poorest;
	}

	public void CapTS(Random flip, double percent, double tax, double alpha, int richpeople, double growth, double gamma)
	{
		 int j=0;
	        int target=0;
	        int richer=0;    // not used, just a label for the richer
	        int poorer=0;
	        double tradeamount=0;
	        double aftertax=0;
	        double totaltax=0;
		    double totalorder=0;
		    growthmatrix=new double[M];
		       
		    sumtrading=0;
		    sumflow=0;
		    int direction=1;  //the parameter to determine the direction of the flow
		    
		    
		    for (int f=0; f< M; f++)
		    {
			   j=flip.nextInt(M);
			   target=findInRange(j, R, flip);    //choose the trading target, then decide who is richer
			   if(wealth[j]>=wealth[target])
			   {
				   richer=j;
				   poorer=target;
				   direction=1;
			   }
			   else
			   {
				   richer=target;
				   poorer=j;
				   direction=-1;
			   }
			   //after deciding who is the poorer, we can decide the trading amount=percent*wealth[poorer]
			   tradeamount=percent*wealth[poorer];
			   aftertax=(1-tax)*tradeamount;
			   totaltax+=tax*tradeamount;
			   
			   
			   sumtrading+=aftertax;
			   
			   //now decide who win this trade
			   if(flip.nextBoolean())    //assume this means j lose
			   {
				   wealth[j]-=aftertax;
				   wealth[target]+=aftertax;
				   sumflow-=direction*aftertax;
			   }
			   else     // this means j wins
			   {
				   wealth[j]+=aftertax;
				   wealth[target]-=aftertax;
				   sumflow+=direction*aftertax;
			   }
			   
		    }
		    
		    totalwealth=0;
		    
		    
		    //now generate the growthmatrix
		    normalization=0;
		    for(int gm=0; gm<M; gm++)
		    {
		    	growthmatrix[gm]=Math.pow(wealth[gm], gamma);
		    	normalization+=growthmatrix[gm];
		    	totalwealth+=wealth[gm];
		    }
		    
		    
		    //now calculate the total growth per step
		    double totalgrowth=0;
		    
		    if(richpeople==M)
		    {
		    	totalgrowth=growth*totalwealth;     //speed up the code if all the agents are considered
		    }
		    else
		    {
		    	double richdata[][]=new double[2][richpeople];
			    
			    richdata=findrich(wealth, richpeople);
			    double totalrichwealth=0;
			    for(int rr=0; rr<richpeople; rr++)
			    {
			    	totalrichwealth+=richdata[1][rr];
			    }
			    
			   totalgrowth= growth*totalrichwealth;
		    }
		    
		    totalwealth=0;   //recalculate the total wealth
		    
			for(int g=0; g<M; g++)   //now everyone's wealth will grow with the amount proportional to growthmatrix[]
			   {
				   wealth[g]=wealth[g]+growthmatrix[g]*totalgrowth/normalization;
				   totalwealth+=wealth[g];  // calculate the total wealth in this loop
			   }
			
		
			   meanwealth=totalwealth/M;
		    
		    
		    for(int ii=0; ii<M; ii++)
		    {
		    	wealth[ii]+=totaltax*(1-alpha)/M;       //after all the trading within one step, everybody gets a benefit from the tax after a dissipation alpha
		    	
		    	if(wealth[ii]!=0)
		    		totalorder-=wealth[ii]/meanwealth*Math.log(wealth[ii]/meanwealth);    //might be problematic if there is tax because of the change in the total wealth
		    }
		    
		    order=totalorder;
	}
	
}

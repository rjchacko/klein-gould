package kang.ising.BasicStructure;


public class Cluster{
	
	public int L1,L2,M;
	public int size;
	public int lattice[];   // the actual distribution of cluster on the lattice
	public int sites[];    // the array of the location labels for the cluster
	public int center; // located the center position of the cluster on the lattice
	
	public Cluster(int L1, int L2, int clustertemp[])
	{
		this.L1=L1;
		this.L2=L2;
		this.M=L1*L2;
		this.lattice= new int [M];
		this.sites= new int [M];
	    this.center=-1;    // -1 means that we have not calculated the position of center yet
	    
		int sizetemp=0;
		for(int i=0; i<M; i++)
		{
			this.sites[i]=clustertemp[i];       // record the locations of the cluster
			this.lattice[i]=0;                  //initialize all the lattice sites
		}
		
		for(int j=0; j<M; j++)
		{
			if(clustertemp[j]>=0)
			{
				this.lattice[clustertemp[j]]=2;
				sizetemp++;
			}
		}
		
		this.size=sizetemp;      //count the number of sites in this cluster
	}
	
	public Cluster()
	{
		this.L1=0;
		this.L2=0;
		this.M=0;
		this.size=0;
		this.center=-1;
		this.sites=null;
		this.lattice=null;
	}
	
	public Cluster clone()
	{
		Cluster copy= new Cluster(L1, L2,sites);
		for(int j=0;j<M;j++)
		{
			copy.lattice[j]=lattice[j];
			copy.sites[j]=sites[j];
		}
		copy.size=size;
		copy.center=center;
		
		return copy;
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
	
	public void Center()    // the function to calculate the center position and restore it in center
	{
		int fp,fx,fy;
		fp=sites[0];
		fx=fp/L2;
		fy=fp%L2;
		double sumx,sumy;
		sumx=0;
		sumy=0;
		int cx,cy;
		int clx[];
		int cly[];
		int tempx,tempy;
		clx= new int[size];
		cly= new int[size];
		clx[0]=fx;
		cly[0]=fy;
		sumx=fx;
		sumy=fy;
		
		for(int cl=1; cl<size; cl++)
		{
			tempx=sites[cl]/L2;
			tempy=sites[cl]%L2;
			clx[cl]=tempx;
			cly[cl]=tempy;
			if((tempx-fx)*(tempx-fx)>(tempx-L1-fx)*(tempx-L1-fx))
				clx[cl]=tempx-L1;
			if((tempx-fx)*(tempx-fx)>(tempx+L1-fx)*(tempx+L1-fx))
				clx[cl]=tempx+L1;
			if((tempy-fy)*(tempy-fy)>(tempy-L2-fy)*(tempy-L2-fy))
				cly[cl]=tempy-L2;
			if((tempy-fy)*(tempy-fy)>(tempy+L2-fy)*(tempy+L2-fy))
				cly[cl]=tempy+L2;
			sumx+=clx[cl];
			sumy+=cly[cl];
		}
		cx=(int)(sumx/size);
		cy=(int)(sumy/size);
		
		
		center=X(cx)*L2+Y(cy);
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
	
	public void CenterDisplay(int display[])
	{
		Center();
		display[center]=3;
		for(int a=0; a<4; a++)
			display[Nneighber(a,center)]=3;
		
	}
 	
 	
}
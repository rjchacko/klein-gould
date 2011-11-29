package kang.ising.BasicStructure;


public class Cluster3D{
	
	public int L1,L2,L3,M;
	public int cx, cy, cz;      //x,y,z indices for display purpose
	
	
	public int size;
	public int lattice[];   // the actual distribution of cluster on the lattice
	public int sites[];    // the array of the location labels for the cluster
	public int center; // located the center position of the cluster on the lattice
	
	public Cluster3D(int L1, int L2, int L3, int clustertemp[])
	{
		this.L1=L1;
		this.L2=L2;
		this.L3=L3;
		this.M=L1*L2*L3;
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
	
	public Cluster3D()
	{
		this.L1=0;
		this.L2=0;
		this.L3=0;
		this.M=0;
		this.size=0;
		this.center=-1;
		this.sites=null;
		this.lattice=null;
	}
	
	public Cluster3D clone()
	{
		Cluster3D copy= new Cluster3D(L1, L2, L3,sites);
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
	
	public int Z(int bz)
	{
		int realz=bz;
		if (bz>=L3)
			realz=bz-L3;
		if (bz<0)
			realz=bz+L3;
		return realz;
	}
	
	public void Center()    // the function to calculate the center position and restore it in center
	{
		int fp,fx,fy, fz;
		fp=sites[0];
		fx=fp/L2;
		fy=fp%L2;
		
		fx=(int)(fp/(L2*L3));
		fy=(int)((int)fp%(L2*L3))/L3;
		fz=(int)((int)fp%(L2*L3))%L3;
		double sumx,sumy,sumz;
		sumx=0;
		sumy=0;
		sumz=0;
		
		int clx[];
		int cly[];
		int clz[];
		int tempx,tempy,tempz;
		clx= new int[size];
		cly= new int[size];
		clz= new int[size];
		clx[0]=fx;
		cly[0]=fy;
		clz[0]=fz;
		sumx=fx;
		sumy=fy;
		sumz=fz;
		
		for(int cl=1; cl<size; cl++)
		{
			
			tempx=(int)(sites[cl]/(L2*L3));
			tempy=(int)((int)sites[cl]%(L2*L3))/L3;
			tempz=(int)((int)sites[cl]%(L2*L3))%L3;
			clx[cl]=tempx;
			cly[cl]=tempy;
			clz[cl]=tempz;
			if((tempx-fx)*(tempx-fx)>(tempx-L1-fx)*(tempx-L1-fx))
				clx[cl]=tempx-L1;
			if((tempx-fx)*(tempx-fx)>(tempx+L1-fx)*(tempx+L1-fx))
				clx[cl]=tempx+L1;
			if((tempy-fy)*(tempy-fy)>(tempy-L2-fy)*(tempy-L2-fy))
				cly[cl]=tempy-L2;
			if((tempy-fy)*(tempy-fy)>(tempy+L2-fy)*(tempy+L2-fy))
				cly[cl]=tempy+L2;
			if((tempz-fz)*(tempz-fz)>(tempz-L3-fz)*(tempz-L3-fz))
				clz[cl]=tempz-L3;
			if((tempz-fz)*(tempz-fz)>(tempy+L3-fz)*(tempz+L3-fz))
				clz[cl]=tempz+L3;
			sumx+=clx[cl];
			sumy+=cly[cl];
			sumz+=clz[cl];
		}
		cx=X((int)(sumx/size));
		cy=Y((int)(sumy/size));
		cz=Z((int)(sumz/size));
		
		
		center=cx*L2*L3+cy*L3+cz;
	}
	
	public int Nneighber(int a,int i ){// function for the index of nearest neighbor
		int nx,ny,nz; //index for neighbor
		int ni=0;
		nx=(int)(i/(L2*L3));
		ny=(int)((int)i%(L2*L3))/L3;
		nz=(int)((int)i%(L2*L3))%L3;
		
		if (a==0) {
			ni=nx*L2*L3+(ny-1)*L3+nz;
			if  (ny==0) {
				ni=nx*L2*L3+(ny+L2-1)*L3+nz;
		}
			
		}//(x,y-1,z) up
		
     	if (a==1){
			ni=(nx+1)*L2*L3+ny*L3+nz;
			if  (nx==L1-1) {
				ni=(nx-L1+1)*L2*L3+ny*L3+nz;
			}
			
		}//(x+1,y,z) right
		
		if (a==2){
			ni=nx*L2*L3+(ny+1)*L3+nz;
			if  (ny==L2-1) {
				ni=nx*L2*L3+(ny-L2+1)*L3+nz;
			}
			
		}//(x,y+1,z) down
		
		if (a==3){
			ni=(nx-1)*L2*L3+ny*L3+nz;
			if  (nx==0) {
				ni=(nx+L1-1)*L2*L3+ny*L3+nz;
			}
		}//(x-1,y,z) left
		
		if (a==4){
			ni=nx*L2*L3+ny*L3+nz+1;
			if  (nz==L3-1) {
				ni=nx*L2*L3+ny*L3+nz-L3+1;
			}
		}//(x,y,z+1) up in z
		
		if (a==5){
			ni=nx*L2*L3+ny*L3+nz-1;
			if  (nz==0) {
				ni=nx*L2*L3+ny*L3+nz+L3-1;
			}
		}//(x,y,z-1) down in z
		
		
		return ni;
		
	}
	
	public void CenterDisplay(int display[])
	{
		Center();
		display[center]=3;
		for(int a=0; a<6; a++)
			display[Nneighber(a,center)]=3;
		
	}
 	
 	
}
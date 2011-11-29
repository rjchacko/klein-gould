package kang.ising.BasicStructure;



public class ClusterSet3D{
	

	public int number;
	public int order[];   // the array to record the location of cluster from the largest to the smallest
	
	public int maximumsize;
	public int maximumpin;
	
	public int minimumsize;
	public int minimumpin;
	public Cluster3D set[];
	
	
	public ClusterSet3D(int number)
	{
		this.number=number;
		this.minimumsize=0;
		this.minimumpin=0;
		this.maximumsize=0;
		this.maximumpin=0;
		this.order=new int[number];
		this.set=new Cluster3D[number];
		for(int j=0;j<number;j++)
			this.set[j]=new Cluster3D();
		
	}
	
	public void findminimum()
	{
		int min=set[0].size;
		int minxy=0;
		for(int j=0; j<number; j++)
		{
			if(set[j].size<=min)
			{
				min=set[j].size;
				minxy=j;
			}
		}
		minimumpin=minxy;
		minimumsize=min;
	}
	
	public void findmaximum()
	{
		int max=set[0].size;
		int maxxy=0;
		for(int j=0; j<number; j++)
		{
			if(set[j].size>=max)
			{
				max=set[j].size;
				maxxy=j;
			}
		}
		maximumpin=maxxy;
		maximumsize=max;
	}
	
	public void AddCluster(Cluster3D temp)
	{
		findminimum();
		if(temp.size>minimumsize)
		{
			set[minimumpin]=temp.clone();
		}
		findminimum();
	}
	
	public void ordering()
	{
		int temp;
		for(int t=0; t<number; t++)
			order[t]=t;
		
		for(int i=0; i<number; i++)	
		{
			for(int j=1; j<number; j++)
		    {
			      if(set[order[j]].size>=set[order[j-1]].size)
			      {
				       temp=order[j];
				       order[j]=order[j-1];
				       order[j-1]=temp;
			      }   
		    }
		}
		maximumsize=set[order[0]].size;
		maximumpin=order[0];
	}
	
	public void CentersDisplay(int map[])
	{
		for(int w=0; w<number; w++)
		{
			set[w].CenterDisplay(map);
		}
	}
	
	public void ClustersDisplay(int map[])     //make all the cluster elements in map[] to be green
	{
		for(int j=0;j<number; j++)
		{
			for(int b=0;b<set[j].M; b++)
			{
				if((map[b]!=2)&set[j].lattice[b]==2)
					map[b]=2;
			}
		}
	}
	
	
}
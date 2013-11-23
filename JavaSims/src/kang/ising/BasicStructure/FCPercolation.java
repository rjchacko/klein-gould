package kang.ising.BasicStructure;

import chris.util.Random;

public class FCPercolation {
	
	public int N;
	public double pb;
	public int clusters[];
	
	public FCPercolation(int n)
	{
		this.N=n;
		this.pb=1;   //default value
		this.clusters=new int [n];
	}
	
	public FCPercolation clone()
	{
		FCPercolation copy=new FCPercolation(N);
		
		copy.SetPb(this.pb);
		
		for (int i=0; i<N; i++)
		{
			copy.clusters[i]=clusters[i];
		}
		return copy;
	}
	
	public void SetPb(double pbN)
	{
		this.pb=pbN/N;
	}
	
	public int PossibleBonds(int n, double pb, Random rand)
	{
		int total=0;
		for(int i=0; i<n; i++)
		{
			if (rand.nextDouble()<=pb)
				total++;
		}
		return total;
	}
	
	public void GenerateClusters(int seed)   //generate clusters and put their size information into array clusters[]
	{
		Random rand=new Random(seed);
		
		int ntemp=0;
		int temp=0;
		int rest=N;
		int index=0;
		int ci=0;  //the increasing clusters index
		
		while(rest>0)
		{
			//initialize the first site in the current cluster
			temp=1;   //temp is the number of sites in the current clusters, not the label of the last site!
			index=0;
			rest--;   //take one from the rest and put it the current cluster
			
			ntemp=PossibleBonds(rest, pb, rand);
			temp+=ntemp;    //connect bonds to the seed site
			rest-=ntemp;    //delete the sites from the rest
		
			//if we haven't checked other sites in the temp cluster, we need to generate cluster from these sites	
			while(index<=(temp-1))
			{
				ntemp=PossibleBonds(rest, pb, rand);
				temp+=ntemp;
				rest=rest-ntemp;
				index++;
			}
			clusters[ci]=temp;
			ci++;
		}
	}
	

}
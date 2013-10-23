package kang.util;

import java.io.File;
import java.io.IOException;
import java.text.DecimalFormat;

import javax.imageio.ImageIO;

import scikit.graphics.dim2.Grid;

public class BasicTools{
	
	public BasicTools()
	{
		
	}
	
 	public double Mean(int data[], int size)
 	{
 		double total=0;
 		double mean=0;
 		for(int q=0; q<size; q++)
 		{
 			total+=data[q];
 		}
 		mean=total/size;
 		return mean;
 	}
 	
 	public double Median(int data[], int size)
 	{
 		int tempdata[]= new int[size];
 		double median;
 		int temp;
 		for(int s=0; s<size; s++)
 		{
 			tempdata[s]=data[s];
 		}
 		for(int x=0; x<size; x++)
 			for(int y=0; y<size-x-1;y++)
 			{
 				if(tempdata[y]>=tempdata[y+1])
 					{
 					temp=tempdata[y];
 					tempdata[y]=tempdata[y+1];
 					tempdata[y+1]=temp;
 					}
 			}
 		median=tempdata[size/2];
 		return median;
 	}
 	
 	
 	public double SD(int data[], int size, double mean)
 	{
 		double totalSD=0;
 		double SD=0;
 		for (int p=0; p<size; p++)
 		{
 			totalSD+=((data[p]-mean)*(data[p]-mean));
 		}
 		SD=Math.sqrt(totalSD/size);
 		
 		return SD;
 	}
 	
 	public double Mean(double data[], int size)
 	{
 		double total=0;
 		double mean=0;
 		for(int q=0; q<size; q++)
 		{
 			total+=data[q];
 		}
 		mean=total/size;
 		return mean;
 	}
 	
 	public double Median(double data[], int size)
 	{
 		double tempdata[]= new double[size];
 		double median;
 		double temp;
 		for(int s=0; s<size; s++)
 		{
 			tempdata[s]=data[s];
 		}
 		for(int x=0; x<size; x++)
 			for(int y=0; y<size-x-1;y++)
 			{
 				if(tempdata[y]>=tempdata[y+1])
 					{
 					temp=tempdata[y];
 					tempdata[y]=tempdata[y+1];
 					tempdata[y+1]=temp;
 					}
 			}
 		median=tempdata[size/2];
 		return median;
 	}
 	
 	public int Max(double data[], int size)
 	{
 		double temp=data[0];
 		int Max=0;
 		for(int s=0; s<size; s++)
 		{
 			if(data[s]>=temp)
 			{
 				Max=s;
 				temp=data[s];
 			}
 		}
 		return Max;
 	}
 	
 	public int Min(double data[], int size)
 	{
 		double temp=data[0];
 		int Min=0;
 		for(int s=0; s<size; s++)
 		{
 			if(data[s]<=temp)
 			{
 				Min=s;
 				temp=data[s];
 			}
 		}
 		return Min;
 	}
 	
 	public double SD(double data[], int size, double mean)
 	{
 		double totalSD=0;
 		double SD=0;
 		for (int p=0; p<size; p++)
 		{
 			totalSD+=((data[p]-mean)*(data[p]-mean));
 		}
 		SD=Math.sqrt(totalSD/size);
 		
 		return SD;
 	}
 	
 	public int Sum(int data[], int size)
 	{
 		int total=0;
 		for(int p=0; p<size; p++)
 		{
 			total+=data[p];
 		}
 		return total;
 	}
	
	public void Picture(Grid grid, int index1, int index2, String path)   //function to capture the grid
	{
		
	    DecimalFormat fmt = new DecimalFormat("0000");
			String SaveAs = path+"pic_"+fmt.format(index1)+"_"+fmt.format(index2)+".png";
		try {
			ImageIO.write(grid.getImage(), "png", new File(SaveAs));
		} catch (IOException e) {
			System.err.println("Error in Writing File" + SaveAs);
		}
		
	}
	
	public double Correlation(double[] map1, double[] map2, int N)  //function to calculate the correlation between two maps
	{
		double C=0;
		double avg1=Mean(map1, N);
		double avg2=Mean(map2, N);
		
		double D1=SD(map1, N, avg1);
		double D2=SD(map2, N, avg2);
	
		double totalc=0;
		for(int i=0; i<N; i++)
		{
			totalc+=((map1[i]-avg1)*(map2[i]-avg2));
		}
		
		C=totalc/(D1*D2*N);
		
		return C;
	}
	
	public double Correlation(double[] map1, int[] map2, int N)  //function to calculate the correlation between two maps
	{
		double C=0;
		double avg1=Mean(map1, N);
		double avg2=Mean(map2, N);
		
		double D1=SD(map1, N, avg1);
		double D2=SD(map2, N, avg2);
	
		double totalc=0;
		for(int i=0; i<N; i++)
		{
			totalc+=((map1[i]-avg1)*(map2[i]-avg2));
		}
		
		C=totalc/(D1*D2*N);
		
		return C;
	}
	
	public double Correlation(int[] map1, int[] map2, int N)  //function to calculate the correlation between two maps
	{
		double C=0;
		double avg1=Mean(map1, N);
		double avg2=Mean(map2, N);
		
		double D1=SD(map1, N, avg1);
		double D2=SD(map2, N, avg2);
	
		double totalc=0;
		for(int i=0; i<N; i++)
		{
			totalc+=((map1[i]-avg1)*(map2[i]-avg2));
		}
		
		C=totalc/(D1*D2*N);
		
		return C;
	}

	public double DCorrelation(double[] map1, double[] map2, int N)  //function to calculate the correlation between two maps except for the dilution
	{
		double C=0;
		double total1=0;
		double total2=0;
		int totalnumber=0;
		
		for(int i=0; i<N; i++)
		{
			if(map1[i]!=0)
			{
				total1+=map1[i];
				total2+=map2[i];
				totalnumber++;
			}
		}
		
		
		
		double avg1=total1/totalnumber;
		double avg2=total2/totalnumber;
		
		double D1=0;
		double D2=0;
		double totalc=0;
		for(int i=0; i<N; i++)
		{
			if(map1[i]!=0)
			{
				D1+=((map1[i]-avg1)*(map1[i]-avg1));
				D2+=((map2[i]-avg2)*(map2[i]-avg2));
				totalc+=((map1[i]-avg1)*(map2[i]-avg2));
			}
			
			
		}
		
		C=totalc/Math.sqrt(D1*D2);
		
		return C;
	}
	
	public int[] BubbleSortIndex(double data[], int index[], boolean ascend)  //return the index of data[] in ascend order(ascend=true) or descend order(ascend=false)
	{
		
		int tempindex=0;

		//bubble sort starts from here	
		for(int out=0; out<data.length-1; out++)
		{
			for(int in=0; in<data.length-1-out; in++)
			{
				if(ascend)  //increasing
				{
					if(data[index[in]]>data[index[in+1]])
					{
						tempindex=index[in];
						index[in]=index[in+1];
						index[in+1]=tempindex;
					}
				}
				else  //decreasing
				{
					if(data[index[in]]<data[index[in+1]])
					{
						tempindex=index[in];
						index[in]=index[in+1];
						index[in+1]=tempindex;
					}
				}
			}
		}
		
		return index;
	}
	
	public int[] IndexToRank(int[] index)   //given the index[], generate the rank of each data point for the original data set
	{
		int[] Rank=new int[index.length];
		for(int i=0; i<index.length; i++)
		{
			Rank[index[i]]=i+1;
		}
		return Rank;
	}
	
	public int turningpoint(double[] data, int Ts, int Te, int dT, boolean decrease)
	{
		double sumdata=0;
		int total=0;
		int time=0;
		
		for(int t=(Ts/dT); t<(Te/dT); t++)
		{
			sumdata+=data[t];
			total++;
		}
	    double avg=sumdata/total;
	    
	    if(decrease)//the data is decreasing in time
	    {
	    	while(data[time/dT]>avg)
	    		{
	    		time+=dT;
	    		}
	    }
	    else   //the data is increasing in time
	    {
	    	while(data[time/dT]<avg)
	    		{
	    		time+=dT;
	    		}
	    }
	    
	    
	    return time;

		
	}
	
	public double[] LorentzData(double[] data, int[] index, boolean sorted)   //the function to generate the Lorentz plot data point for data[]
	{
		int N=data.length;
		double Ldata[]=new double [N];
		
		double total=0;   //the sum of data
		
		for(int i=0; i<N; i++)
		{
			total+=data[i];
		}
					
		
		if(sorted)
		{

			// if the original data is already sorted in ascending order, then index[i]=i
		}
		else
		{
		   index=BubbleSortIndex(data, index, true);
		}
		
		//after sorting, calculate the Lorentz data
		
		for(int j=0; j<N; j++)
		{
			Ldata[j]=0;
			
			for(int k=0; k<j; k++)
			{
				Ldata[j]+=(data[index[k]]/total);
			}
				
		}
		
		return Ldata;
	}
	
	public double Gini(double[] Ldata)  //calculate the Gini coefficient from the Lorentz data point
	{
		int N=Ldata.length;
		double G=0;
		double totalLdata=0;
		for(int i=0; i<N; i++)
		{
			totalLdata+=Ldata[i];
		}
		G=1-2*totalLdata/(N+1);
		return G;
	}
}
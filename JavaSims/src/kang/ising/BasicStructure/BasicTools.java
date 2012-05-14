package kang.ising.BasicStructure;

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
	
	
	
	
}
//package kang.ising;
//
//import java.awt.Color;
//
//import chris.util.PrintUtil;
//
//import scikit.graphics.ColorPalette;
//import scikit.graphics.dim2.Grid;
//import scikit.jobs.Control;
//import scikit.jobs.Job;
//import scikit.jobs.Simulation;
//import scikit.jobs.params.DoubleValue;
//
//public class nucleationising extends Simulation{
//	Grid grid1 = new Grid("Diluted Ising lattice 2d for nucleation");
//
//	public int L1, L2, M; //parameters for the lattice
//    public int i, x, y;  //parameters for the index
//	public double T;     //temperature 
//	public double H;     //field
//	public double J;     //interaction constant after normalization
//	public double NJ;    //interaction constant before normalization
//	public double totalmagnetization;
//	public double magnetization;
//	public int R;   //interaction range
//	
//	
//	public int deadsites;
//	public int step; //Monte Carlo step
//	public int steplimit; //upperlimit of MC step
//
//	public int isingspin[];     //the array of the data
//	
//	
//	public int Nneighber(int a,int i ){// function for the index of nearest neighbor
//		int nx,ny; //index for neighbor
//		int ni=0;
//		nx=(int)i/L2;
//		ny=(int)i%L2;
//		
//		if (a==0) {
//			ni=nx*L2+ny-1;
//			if  (ny==0) {
//				ni=nx*L2+ny+L2-1;
//			}
//			
//		}//(x,y-1) up
//		
//		if (a==1){
//			ni=(nx+1)*L2+ny;
//			if  (nx==L1-1) {
//				ni=(nx-L1+1)*L2+ny;
//			}
//			
//		}//(x+1,y) right
//		
//		if (a==2){
//			ni=nx*L2+ny+1;
//			if  (ny==L2-1) {
//				ni=nx*L2+ny-L2+1;
//			}
//			
//		}//(x,y+1) down
//		
//		if (a==3){
//			ni=(nx-1)*L2+ny;
//			if  (nx==0) {
//				ni=(nx+L1-1)*L2+ny;
//			}
//		}//(x-1,y) left
//		
//		return ni;
//		
//	}
//	
//	public double interactionE (int j){ //function for interaction energy
//		double Energy=0;
//		int b,k;
//		for(b=0; b<4;b++){
//			k=Nneighber(b,j);
//			Energy=Energy+J*isingspin[j]*isingspin[k];
//		}
//		return Energy;
//			
//	}
//	
//	public double longrangeE (int i){
//		double Energy=0;
//		int S=0;
//		int nx=i/L2;
//		int ny=i%L2;
//		int kx, ky;
//		
//		for (int m=-R; m<=R; m++)
//			for (int n=-R; n<=R; n++)
//			{
//				kx=nx+m;
//				ky=ny+n;
//				if(nx+m<0)
//					kx=nx+m+L1;
//				if(nx+m>L1-1)
//					kx=nx+m-L1;
//				if(ny+n<0)
//					ky=ny+n+L2;
//				if(ny+n>L2-1)
//					ky=ny+n-L2;
//				S+=isingspin[kx*L2+ky];	
//			}
//		Energy=J*isingspin[i]*S-J;
//		return Energy;
//	}
//	
//	public static void main (String[] kangliu){
//		new Control(new nucleationising(), "Kang Liu's diulted ising model for nucleation" );
//	}
//	
//	
//	
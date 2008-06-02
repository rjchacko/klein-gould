package chris.ofc;


import java.io.BufferedInputStream;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.text.DecimalFormat;

import javax.imageio.ImageIO;

import scikit.graphics.dim2.Grid;
import scikit.jobs.Job;
import scikit.jobs.params.Parameters;
import chris.util.DirUtil;
import chris.util.LatticeNeighbors;
import chris.util.PrintUtil;

public class NfailDamage2D extends SimpleDamage2D{
	
	// Parameters
	 public double Sr0, Sc0, Srwidth, Scwidth, alphawidth, lifewidth, DaboutFS, tkip, dtkip;
	 public int Nlives, rmin, hammersize, Nshowers, SonFSindex, NdeadS, Nalive, search;
	 public String lifeshape, residualnoise, criticalnoise, outdir, outfile1, PicDir, outfile2;
	 public double Sr[], Sc[], SsoFar[];
	 public Boolean ShowGrid;
	 public double SonFS[], NlivesLeft[];
	 
	 // Formats
	 
	 public DecimalFormat fmt = new DecimalFormat("0000000");
	 public DecimalFormat fmts = new DecimalFormat("0000");

	 //String outfileCHK;
	 
	// Constructor
	public NfailDamage2D(Parameters params) {
		
		super(params); // I assume this is a call to SimpleDamage2D's constructor
	
		PseudoConstructorNF(params);
		
	}
	
	
	public void PseudoConstructorNF(Parameters params){
	
		Sr0           = params.fget("Residual Stress (\u03C3_r)");
		Srwidth       = params.fget("\u03C3_r width");
		Sc0           = params.fget("Critical Stress (\u03C3_c)");
		Scwidth       = params.fget("\u03C3_c width");
		lifewidth     = params.fget("Nlives Width");
		Nlives        = params.iget("Number of Lives");
		rmin          = params.iget("Minimum Interaction Radius (r)");	
		lifeshape     = params.sget("Life Style");
		residualnoise = params.sget("\u03C3_r Noise");
		criticalnoise = params.sget("\u03C3_r Noise");
		outdir        = params.sget("Data Directory");

		outfile1=outdir+File.separator+"Damage1.txt";
		outfile2=outdir+File.separator+"Damage2.txt";
		//outfileCHK=outdir+File.separator+"DamageCHK.txt";
		
		PicDir=outdir+"/Pics/";
		DirUtil.MkDir(PicDir);
		
		alphawidth = Nwidth;
		
		alive      = new int[2*N];
		Sr         = new double[N];
		Sc         = new double[N];
		SsoFar     = new double[N];
		SonFS      = new double[Nlives*N];
		NlivesLeft = new double[Nlives+1];
		
		if (Nlives == 1){
			Sr0=0.;
			Srwidth=0.;
		}
		
		return;
	}
	
	public void Initialize(String str){
		
		
		
		if (BCs.equals("Bordered")){
			if(shape.equals("Circle")){
				neighbors = new LatticeNeighbors(L,L,rmin,R,LatticeNeighbors.Type.BORDERED,LatticeNeighbors.Shape.Circle);
			}
			else if(shape.equals("Square")){
				neighbors = new LatticeNeighbors(L,L,rmin,R,LatticeNeighbors.Type.BORDERED,LatticeNeighbors.Shape.Square);
			}
			else if(shape.equals("Diamond")){
				neighbors = new LatticeNeighbors(L,L,rmin,R,LatticeNeighbors.Type.BORDERED,LatticeNeighbors.Shape.Diamond);
			}
		}
		else{
			if(shape.equals("Circle")){
				neighbors = new LatticeNeighbors(L,L,rmin,R,LatticeNeighbors.Type.PERIODIC,LatticeNeighbors.Shape.Circle);
			}
			else if(shape.equals("Square")){
				neighbors = new LatticeNeighbors(L,L,rmin,R,LatticeNeighbors.Type.PERIODIC,LatticeNeighbors.Shape.Square);
			}
			else if(shape.equals("Diamond")){
				neighbors = new LatticeNeighbors(L,L,rmin,R,LatticeNeighbors.Type.PERIODIC,LatticeNeighbors.Shape.Diamond);
			}
		}
				
		if(str.equals("Flat")){
			
			for (int i = 0 ; i < N ; i++){

				if(criticalnoise.equals("On")) {
					Sc[i]=Scwidth*rand.nextGaussian()+Sc0;
				}
				else{
					Sc[i]=Sc0;
				}
				
				stress[i]  = Sc0*rand.nextDouble();
				if((Sc[i]-stress[i])<(Sc[imax]-stress[imax])) imax=i;
				alive[i+N] = 1;
				
				if (residualnoise.equals("On")) {
					Sr[i]=Srwidth*rand.nextGaussian()+Sr0;
				}
				else{
					Sr[i]=Sr0;
				}
				
				if (lifeshape.equals("Flat")){
					alive[i]=rand.nextInt(Nlives)+1;
				}
				else if (lifeshape.equals("Gaussian")){
					alive[i]=(int)(Math.floor((lifewidth*rand.nextGaussian()+1.*Nlives)+0.5));
				}
				else if(lifeshape.equals("Constant")){
					alive[i]=Nlives;
				}
			}
			
			for (int ii = 0 ; ii < Nlives ; ii++){
				NlivesLeft[ii] = 0;
			}
			NlivesLeft[Nlives] = N;
			
			NlivesLeft[alive[imax]]--;
			NlivesLeft[alive[imax]-1]++;
			
			alive[imax]--;
			alive[imax+N]=0;
			
			// Bring Site with Most Stress to Failure
			
			stressMax = stress[imax];
			for (int i = 0; i<N; i++){
				stress[i]+=Sc[imax]-stressMax;
				SsoFar[i]=stress[i];
			}
			dtkip = Sc[imax] - stressMax;
			
		}
		else {
			System.out.println("Error! Intialization type " + str + " does not exist!");
		}
		
		time=0;
		tkip=0;
		showernumber=0;
				
		return;
	}	
	
	public void Initialize(Parameters prms){
		
		// Set everything up as if the method Avalanche() just returned
		
		// Find the directory in which the data files live
		String InDir = prms.sget("Input Directory");
		
		// Read in the parameters
		String ParamsIn = InDir + File.separator + "Params4Clone.txt";
		
		try{
			String record;
			int posDelimiter;

			File f = new File(ParamsIn); 
			FileInputStream fis = new FileInputStream(f); 
			BufferedInputStream bis = new BufferedInputStream(fis); 
			BufferedReader bir = new BufferedReader(new InputStreamReader(bis));	

			// Make output directory
			record = bir.readLine();
			posDelimiter = record.indexOf('='); 
			String OD = record.substring(posDelimiter + 2);
			OD += "/Continue";
			DirUtil.MkDir(OD);
			prms.set("Data Directory",OD);

			// Read in the rest of the parameters and set them
			record = bir.readLine();
			posDelimiter = record.indexOf('='); 
			prms.set("Random Seed",(int)(Double.parseDouble(record.substring(posDelimiter + 2))));

			record = bir.readLine();
			posDelimiter = record.indexOf('='); 
			prms.set("Animation",record.substring(posDelimiter + 2));

			record = bir.readLine();
			posDelimiter = record.indexOf('='); 
			prms.set("Lattice Size",(int)(Double.parseDouble(record.substring(posDelimiter + 2))));

			record = bir.readLine();
			posDelimiter = record.indexOf('='); 
			prms.set("Number of Lives",(int)(Double.parseDouble(record.substring(posDelimiter + 2))));

			record = bir.readLine();
			posDelimiter = record.indexOf('='); 
			prms.set("Life Style",record.substring(posDelimiter + 2));

			record = bir.readLine();
			posDelimiter = record.indexOf('='); 
			prms.set("Nlives Width",(Double.parseDouble(record.substring(posDelimiter + 2))));

			// skip T_max
			record = bir.readLine();

			record = bir.readLine();
			posDelimiter = record.indexOf('=');  
			prms.set("Boundary Condtions",record.substring(posDelimiter + 2));

			record = bir.readLine();
			posDelimiter = record.indexOf('='); 
			prms.set("Critical Stress (\u03C3_c)",(Double.parseDouble(record.substring(posDelimiter + 2))));

			record = bir.readLine();
			posDelimiter = record.indexOf('='); 		
			prms.set("\u03C3_c Noise",record.substring(posDelimiter + 2));

			record = bir.readLine();
			posDelimiter = record.indexOf('='); 		
			prms.set("\u03C3_c width",(Double.parseDouble(record.substring(posDelimiter + 2))));

			record = bir.readLine();
			posDelimiter = record.indexOf('='); 		
			prms.set("Residual Stress (\u03C3_r)",(Double.parseDouble(record.substring(posDelimiter + 2))));

			record = bir.readLine();
			posDelimiter = record.indexOf('='); 		
			prms.set("\u03C3_r Noise", record.substring(posDelimiter + 2));

			record = bir.readLine();
			posDelimiter = record.indexOf('='); 		
			prms.set("\u03C3_r width",(Double.parseDouble(record.substring(posDelimiter + 2))));

			record = bir.readLine();
			posDelimiter = record.indexOf('='); 	
			prms.set("Interaction Shape", record.substring(posDelimiter + 2));

			record = bir.readLine();
			posDelimiter = record.indexOf('='); 		
			prms.set("Interaction Radius (R)",(int)(Double.parseDouble(record.substring(posDelimiter + 2))));

			record = bir.readLine();
			posDelimiter = record.indexOf('='); 		
			prms.set("Minimum Interaction Radius (r)",(int)(Double.parseDouble(record.substring(posDelimiter + 2))));

			record = bir.readLine();
			posDelimiter = record.indexOf('='); 	
			prms.set("Dissipation (\u03B1)",(Double.parseDouble(record.substring(posDelimiter + 2))));

			record = bir.readLine();
			posDelimiter = record.indexOf('='); 		
			prms.set("\u03B1 Noise",record.substring(posDelimiter + 2));

			record = bir.readLine();
			posDelimiter = record.indexOf('='); 	
			prms.set("\u03B1 Width",(Double.parseDouble(record.substring(posDelimiter + 2))));

			record = bir.readLine();
			posDelimiter = record.indexOf('='); 	
			prms.set("Record",record.substring(posDelimiter + 2));

			record = bir.readLine();
			posDelimiter = record.indexOf('='); 	
			prms.set("Number of Resets",(Double.parseDouble(record.substring(posDelimiter + 2))));

			record = bir.readLine();
			posDelimiter = record.indexOf('='); 	
			prms.set("Number of Showers",(int)(Double.parseDouble(record.substring(posDelimiter + 2))));
			
			// Print new params
			PrintParams(OD + File.separator + "Params4Continue.txt", prms);
			
			// Reset all parameters that are set by the constructor(s)
		
			PseudoConstructorSD(prms);
			PseudoConstructorNF(prms);
			
			// Set Up Lattice Neighbors Shape
			
			if (BCs.equals("Bordered")){
				if(shape.equals("Circle")){
					neighbors = new LatticeNeighbors(L,L,rmin,R,LatticeNeighbors.Type.BORDERED,LatticeNeighbors.Shape.Circle);
				}
				else if(shape.equals("Square")){
					neighbors = new LatticeNeighbors(L,L,rmin,R,LatticeNeighbors.Type.BORDERED,LatticeNeighbors.Shape.Square);
				}
				else if(shape.equals("Diamond")){
					neighbors = new LatticeNeighbors(L,L,rmin,R,LatticeNeighbors.Type.BORDERED,LatticeNeighbors.Shape.Diamond);
				}
			}
			else{
				if(shape.equals("Circle")){
					neighbors = new LatticeNeighbors(L,L,rmin,R,LatticeNeighbors.Type.PERIODIC,LatticeNeighbors.Shape.Circle);
				}
				else if(shape.equals("Square")){
					neighbors = new LatticeNeighbors(L,L,rmin,R,LatticeNeighbors.Type.PERIODIC,LatticeNeighbors.Shape.Square);
				}
				else if(shape.equals("Diamond")){
					neighbors = new LatticeNeighbors(L,L,rmin,R,LatticeNeighbors.Type.PERIODIC,LatticeNeighbors.Shape.Diamond);
				}
			}

			// Read in variable value (e.g. stess[])
			
			// Need five data readers
			
			String rr1, rr2, rr3, rr4, rr5;
			int pd1, pd2, pd3, pd4, pd5;

			File f1 = new File(InDir + File.separator + "Alive.txt"); 
			File f2 = new File(InDir + File.separator + "Critical_Stress.txt"); 
			File f3 = new File(InDir + File.separator + "Residual_Stress.txt"); 
			File f4 = new File(InDir + File.separator + "Stress_so_Far.txt"); 
			File f5 = new File(InDir + File.separator + "Stress.txt"); 
			
			FileInputStream fis1 = new FileInputStream(f1); 
			BufferedInputStream bis1 = new BufferedInputStream(fis1); 
			FileInputStream fis2 = new FileInputStream(f2); 
			BufferedInputStream bis2 = new BufferedInputStream(fis2); 
			FileInputStream fis3 = new FileInputStream(f3); 
			BufferedInputStream bis3 = new BufferedInputStream(fis3); 
			FileInputStream fis4 = new FileInputStream(f4); 
			BufferedInputStream bis4 = new BufferedInputStream(fis4); 
			FileInputStream fis5 = new FileInputStream(f5); 
			BufferedInputStream bis5 = new BufferedInputStream(fis5); 
						
			BufferedReader b1 = new BufferedReader(new InputStreamReader(bis1));	
			BufferedReader b2 = new BufferedReader(new InputStreamReader(bis2));
			BufferedReader b3 = new BufferedReader(new InputStreamReader(bis3));
			BufferedReader b4 = new BufferedReader(new InputStreamReader(bis4));
			BufferedReader b5 = new BufferedReader(new InputStreamReader(bis5));

			// skip first ten (10) lines
			
			for (int headerskip = 0 ; headerskip < 10 ; headerskip++){
			
				rr1 = b1.readLine();
				rr2 = b2.readLine();
				rr3 = b3.readLine();
				rr4 = b4.readLine();
				rr5 = b5.readLine();
			
			}
			
			// start reading in values
			
			for (int jj = 0 ; jj < L ; jj++){
				
				rr1 = b1.readLine();
				rr2 = b2.readLine();
				rr3 = b3.readLine();
				rr4 = b4.readLine();
				rr5 = b5.readLine();			
				pd1 = rr1.indexOf('\t'); 
				pd2 = rr2.indexOf('\t'); 
				pd3 = rr3.indexOf('\t'); 
				pd4 = rr4.indexOf('\t'); 
				pd5 = rr5.indexOf('\t'); 
				
				for (int kk = 0 ; kk < L-1 ; kk++){
					
					alive[jj*L+kk]  = (int)(Double.parseDouble(rr1.substring(0, pd1)));
					Sc[jj*L+kk]     =       Double.parseDouble(rr2.substring(0, pd2));
					Sr[jj*L+kk]     =       Double.parseDouble(rr3.substring(0, pd3));
					SsoFar[jj*L+kk] =       Double.parseDouble(rr4.substring(0, pd4));
					stress[jj*L+kk] =       Double.parseDouble(rr5.substring(0, pd5));
					
					rr1 = rr1.substring(pd1 + 1);
					rr2 = rr2.substring(pd2 + 1);
					rr3 = rr3.substring(pd3 + 1);
					rr4 = rr4.substring(pd4 + 1);
					rr5 = rr5.substring(pd5 + 1);
					pd1 = rr1.indexOf('\t'); 
					pd2 = rr2.indexOf('\t'); 
					pd3 = rr3.indexOf('\t'); 
					pd4 = rr4.indexOf('\t'); 
					pd5 = rr5.indexOf('\t');

				}
				
				// read in last value in row jj
				
				alive[jj*L+L-1]  = (int)(Double.parseDouble(rr1));
				Sc[jj*L+L-1]     =       Double.parseDouble(rr2);
				Sr[jj*L+L-1]     =       Double.parseDouble(rr3);
				SsoFar[jj*L+L-1] =       Double.parseDouble(rr4);
				stress[jj*L+L-1] =       Double.parseDouble(rr5);
				
			}
			
			// continue on for alive[] only
			
			for (int jj = 0 ; jj < L ; jj++){
				
				rr1 = b1.readLine();
				pd1 = rr1.indexOf('\t'); 
				
				for (int kk = 0 ; kk < L - 1; kk++){
					
					alive[jj*L+kk+N]  = (int)(Double.parseDouble(rr1.substring(0, pd1)));
					rr1 = rr1.substring(pd1 + 1);
					pd1 = rr1.indexOf('\t'); 
			
				}
				
				alive[jj*L+L-1+N]  = (int)(Double.parseDouble(rr1));
				
			}
			
			// read in last line from data file
			// set things like time etc.
			
			File fL = new File(InDir + File.separator + "Data_LL.txt");
			FileInputStream fisL = new FileInputStream(fL); 
			BufferedInputStream bisL = new BufferedInputStream(fisL); 
			BufferedReader bL = new BufferedReader(new InputStreamReader(bisL));
			
			String rrL;
			int pdL;
			
			// skip first line
			rrL = bL.readLine();			
			pdL = rrL.indexOf('\t');
			
			rrL = bL.readLine();			
			pdL = rrL.indexOf('\t'); 
			
			time = Double.parseDouble(rrL.substring(0, pdL));
			rrL = rrL.substring(pdL + 1);
			tkip = Double.parseDouble(rrL.substring(0, pdL));
			rrL = rrL.substring(pdL + 1);
			Nshowers = (int) Double.parseDouble(rrL.substring(0, pdL));
			rrL = rrL.substring(pdL + 1);
			NdeadS = (int) Double.parseDouble(rrL.substring(0, pdL));
			
			showernumber=0;

			// set imax
			
			imax = 0;
			for (int kk = 0 ; kk < N ; kk++){
				if(stress[kk] > stress[imax]) imax = kk;
			}
			
			
		}
		catch (IOException e) { 
			System.out.println("ERROR!" + e.getMessage()); 
			crack=true; // kill app
		}
		
		return;
	}

	public void Avalanche() {
		
		SonFSindex=0;
		showernumber=0;
		Nshowers=0;
		
		time++;
		tkip+=dtkip;
		
		if(!(shape.equals("All Sites"))){
			DaboutFS = FSdensity(imax);
		}
		else{
			DaboutFS = -1;
		}
		
		// Distribute stress from initial failure
		
		dead[0]=imax;
		search = 1;
		DistStress(shape);
		
		// reset plate
		resetPlate();
		
		// search for an avalanche
		findAvalanche();
		
		while (search>0){
			
			showernumber++;
			// redistribute stress of failure(s)
			DistStress(shape);
			
			// reset plates
			resetPlate();
			
			// look for subsequent avalanche
			findAvalanche();
			
			showering=true;
			Job.animate();
	
		}
		// set up next failure

		
		// Find most stressed site
		imax=0;
		for (int i = 0 ; i < N ; i++){
			if((Sc[i]-stress[i])<(Sc[imax]-stress[imax])) imax=i;
		}
		
		// check if most stressed site has already failed
		if(alive[imax] == 0){
			System.out.println("All Sites Failed!");
			crack=true;
			return;
		}
		
		// kill the site 
		NlivesLeft[alive[imax]]--;
		NlivesLeft[alive[imax]-1]++;
		alive[imax]--;
		alive[imax+N]=0;
		
		// bring the most stressed site to failure
		stressMax = stress[imax];
		for (int i = 0; i<N; i++){
			stress[i]+=Sc[imax]-stressMax;
		}
		dtkip = Sc[imax] - stressMax;
		
		showering=false;
		Job.animate();
				
		return;
	}

	public void DistStress(String HowToDumpStress){		
		
		if(HowToDumpStress.equals("All Sites")){
			
			for (int i = 0 ; i < search ; i++){
				Nalive = 0;
				for (int j = 0; j<N; j++){
					Nalive+=alive[j+N];
				}	
				if (Nalive > 0){
					if(alive[dead[i]]>0){
						release=(stress[dead[i]]-Sr[dead[i]])/Nalive;
					}
					else{
						release=stress[dead[i]]/Nalive;
					}
					
					if(Nbool){
						for (int j = 0 ; j < dead[i] ; j++){
							stress[j]+=(1-alphawidth*rand.nextGaussian()-alpha)*release*alive[j+N];
						}
						for (int j = dead[i]+1 ; j < N ; j++){
							stress[j]+=(1-alphawidth*rand.nextGaussian()-alpha)*release*alive[j+N];
						}
					}
					else{
						for (int j = 0 ; j < dead[i] ; j++){
							stress[j]+=(1-alpha)*release*alive[j+N];
						}
						for (int j = dead[i]+1 ; j < N ; j++){
							stress[j]+=(1-alpha)*release*alive[j+N];
						}
					}
				}
			}
		}
		else{	// distribute the stress to a subset of the sites on the lattice
			
			for (int i = 0; i<search; i++){

				int[] nbs = neighbors.get(dead[i]);
				Nalive=0;
				for (int j = 0; j<nbs.length; j++){
					Nalive+=alive[nbs[j]+N];
				}			
				if(Nalive>0){									
					if(alive[dead[i]]>0){
						release=(stress[dead[i]]-Sr[dead[i]])/Nalive;
					}
					else{
						release=stress[dead[i]]/Nalive;
					}
					
					
					if(Nbool){
						for (int j = 0; j<nbs.length; j++){
							stress[nbs[j]]+=(1-alphawidth*rand.nextGaussian()-alpha)*release*alive[nbs[j]+N];
						}
					}
					else{
						for (int j = 0; j<nbs.length; j++){
							stress[nbs[j]]+=(1-alpha)*release*alive[nbs[j]+N];
						}
					}
				}
			}
		}
		
		return;
	}
	
	public int[] getStressLines(int center, int Nlines){
		
		// NB this method *assumes* a circular potential
		
		int[] temp = new int[Nlines*(2*R+2)];
		int[] ret;
		int count = 0;
		
		// Generate the trig fns
		
		double sintheta = 2*rand.nextDouble() - 1;
		double costheta = Math.sqrt(1-sintheta*sintheta);
		double tantheta = sintheta/costheta;
		
		if(BCs.equals("Bordered")){
		
			int x0 = center%L;
			int y0 = (int)(center/L);
			int minG, maxG, xg, yg;
			
			for (int countlines = 0 ; countlines < Nlines ; countlines++){

				// Separate cases for  | slope | > / < 1
				
				if(costheta > 1/Math.sqrt(2)){		// |slope| < 1
					
					minG = x0 - (int)(R*costheta);
					if (minG < 0) minG = 0;
		
					maxG = x0 + (int)(R*costheta);
					if (maxG > L) maxG = L;
				
					for (int ii = 0 ; ii < maxG - minG + 1 ; ii++){
						xg = minG + ii;
						yg = (int)(Math.round(tantheta*(xg-x0) + y0));
						
						if (yg > 0 && yg < L) temp[count++] = yg*L + xg;
					}
				}
				else{	// |slope| > 1
					
					minG = y0 - Math.abs((int)(R*sintheta));
					if (minG < 0) minG = 0;
					
					maxG = y0 + Math.abs((int)(R*sintheta));
					if (maxG > L) maxG = L;
					
					for (int ii = 0 ; ii < maxG - minG + 1 ; ii++){
						yg = minG + ii;
						xg = (int)(Math.round((yg-y0)/tantheta + x0));
						
						if (xg > 0 && xg < L) temp[count++] = yg*L + xg;
					}	
				}
			}
		}
		else{		// BCs are periodic
			
			int x0 = center%L;
			int y0 = (int)(center/L);
			int minl, xl, yl;
			
			for (int countlines = 0 ; countlines < Nlines ; countlines++){
				
				if(costheta > 1/Math.sqrt(2)){		// |slope| < 1
				
					minl = x0 - (int)(R*costheta);
					
					for (int ii = 0 ; ii < (int)(2*R*costheta)+1 ; ii++){
						xl = minl + ii;
						yl = (int)(Math.round(tantheta*(xl-x0) + y0));
						
						if(xl < 0){
							xl+=L;
						}
						else if(xl > L){
							xl = L - xl;
						}
						if(yl < 0){
							yl+=L;
						}
						else if(yl > L){
							yl = L - yl;
						}
						
						temp[count++] = yl*L + xl;
						
					}
					
					
				}
				else{	// |slope| > 1
					
					minl = y0 - Math.abs((int)(R*sintheta));
					
					for (int ii = 0 ; ii < (int)(2*Math.abs((int)(R*sintheta)))+1 ; ii++){
						yl = minl + ii;
						xl = (int)(Math.round((yl-y0)/tantheta + x0));
						
						if(xl < 0){
							xl+=L;
						}
						else if(xl > L){
							xl = L - xl;
						}
						if(yl < 0){
							yl+=L;
						}
						else if(yl > L){
							yl = L - yl;
						}
						
						temp[count++] = yl*L + xl;
						
					}
				}
			}
		}
		
		ret = new int[count];
		
		for (int ii = 0 ; ii < count ; ii++){
			ret[ii] = temp[ii];
		}
		
		
		return ret;
	}
	
	public void resetPlate(){

		for (int i = 0; i<search; i++){
			if (alive[dead[i]]>0) {
				stress[dead[i]]=Sr[dead[i]];
				alive[dead[i]+N]=1;
			}
			else {
				stress[dead[i]]=-2;
			}
		}
		
		return;
	}
	
	public void findAvalanche(){
		
		// can really reduce the search to the neighbor sites
		
		search=0;
		
		for (int i = 0; i<N ; i++){
			if(stress[i]>Sc[i]){
				SonFS[SonFSindex++]=stress[i];
				dead[search++]=i;
				alive[i+N]=0;
				NlivesLeft[alive[i]]--;
				if(alive[i]>0) NlivesLeft[alive[i]-1]++;
				alive[i]--;
			}
		}
		
		Nshowers+=search;
		
		return;
	}

	public int[] LiveSites(){
		
		NdeadS = 0;
		
		int[] tmp, rtn;
		int j=0, nonzero=0;
		
		tmp = new int[N];
		
		for( int i = 0 ; i < N ; i++) {
			tmp[i] = alive[i+N];
			if (alive[i+N] != 0) nonzero++;
			if (alive[i] == 0) NdeadS++;
		}
		
		rtn = new int[nonzero];
		
		for ( int i = 0 ; i < N ; i++){
			if(tmp[i] != 0){
				rtn[j++]=i;
			}
		}
		
		return rtn;
	}
	
	public double GetMax(double[] array){
		
		double max=array[0];
		
		for (int i = 0 ; i < array.length ; i++){
			if(array[i]>max) max=array[i];
		}
		
		
		return max;
	}
	
	public int GetMax(int[] array){
		
		int max=array[0];
		
		for (int i = 0 ; i < array.length ; i++){
			if(array[i]>max) max=array[i];
		}
		
		
		return max;
	}
	
	public double radiusGyration(int index){
		
		if(shape.equals("All Sites")) return -1;
		
		long x0,y0,x,y,dx,dy;
		long mass=0;
		double rg=0;
		
		int[] RGnbs = neighbors.get(index);
		
		x0=index%L;
		y0=(int)(index/L);
		if (BCs.equals("Bordered")){
			for (int i=0 ; i<RGnbs.length ; i++){
				x=RGnbs[i]%L;
				y=(int)(RGnbs[i]/L);
				mass+=alive[RGnbs[i]+N];
				dx=(x-x0);
				dy=(y-y0);
				rg+=alive[RGnbs[i]+N]*(dx*dx+dy*dy);
			}
		}
		else{
			for (int i=0 ; i<RGnbs.length ; i++){
				x=RGnbs[i]%L;
				y=(int)(RGnbs[i]/L);
				mass+=alive[RGnbs[i]+N];
				dx=(x-x0);
				dy=(y-y0);
				if (2*Math.abs(dx) > L) dx=L-Math.abs(dx);
				if (2*Math.abs(dy) > L) dy=L-Math.abs(dy);
				rg+=alive[RGnbs[i]+N]*(dx*dx+dy*dy);
			}

		}
		
		if(mass>0){
			rg=rg/mass;
		}
		rg=Math.sqrt(rg);
		
		return rg;
	}
	
	public double EFmetric(){
		
		double ret=0;	
		double Sbar=0;
		
		
		for (int i = 0 ; i < N ; i++){
			SsoFar[i]+=stress[i];
		}
		
		for (int i = 0 ; i < N ; i++){
			Sbar+=SsoFar[i]*alive[i+N];
		}
		
		Sbar=Sbar/N;
		
		for (int i = 0 ; i < N ; i++){
			ret+=(SsoFar[i]*alive[i+N]-Sbar)*(SsoFar[i]*alive[i+N]-Sbar);
		}
		
		ret=ret/(time*time);
		
		return ret;
	}
	
	public double FSdensity(int index){
		
		double ret=0;
		int[] DSnbs = neighbors.get(index);
		
		if (DSnbs.length > 1){
			
			for (int i = 0 ; i < DSnbs.length ; i++){
				ret+=alive[DSnbs[i]+N];
			}
			
			// remove failure of imax 
			ret--;
			
			// divide by all sites save imax
			ret=ret/(DSnbs.length-1);
			
			// ret is the density of life site so use 
			// rho_live + rho_dead = 1
			
			ret=1-ret;

		}
		return ret;
	}
	
	public double GetAve(double set[]){

		int temp = set.length;
	
		return GetAve(set,temp);
	}
	
	public double GetAve(double set[], int Lset){

		double ret=0;
		
		if (Lset > 0){
		
			for (int i = 0 ; i < Lset ; i++){
				ret+=set[i];
			}
			
			ret=ret/Lset;
		
		}
		
		return ret;
	}
	
	public void WriteDataHeader(){
	
		WriteDataHeader(outfile1,1);
		WriteDataHeader(outfile2,2);
		
		return;
	}
	
	public void WriteDataHeader(String fout, int Num){
		
		switch(Num){
		
			case 1:
				PrintUtil.printlnToFile(fout,"Time","t_kip","N_avlnchs","N_dead","Rgyr","Omega","<FS_stress>","rho_FS");
				break;
			case 2:
				PrintUtil.printlnToFile(fout,"Time","Nlives=0","Nlives=1","Nlives=2",". . .","Nlives=Nmax");
				break;
			default:
				System.err.println("File Not Found!");
				break;
		}
		return;
	}
	
	public void TakeData(){
		
		TakeData(outfile1,outfile2);
		
		return;
	}
	
	public void TakeData(String fout1, String fout2){
		
		double rgyr;

		int[] LS = LiveSites(); 
		int LSlength = LS.length;
		if(LSlength>0){
			rgyr=radiusGyration(LS[rand.nextInt(LSlength)]);
		}
		else{
			rgyr=0;
		}

		PrintUtil.printlnToFile(fout1,time,tkip,Nshowers,NdeadS,rgyr,EFmetric(),GetAve(SonFS,SonFSindex),DaboutFS);
		PrintUtil.printTimeAndVectorToFile(fout2, time, NlivesLeft);


		return;
	}
	
	public static void PrintParams(String fout, Parameters prms){
		
		PrintUtil.printlnToFile(fout,prms.toString());
		
		return;
	}
	
	public void TakePicture(Grid grid){
		
		String SaveAs = PicDir + File.separator + grid.getTitle()+fmt.format(time)+"-"+fmts.format(showernumber)+".png";
		
		try {
			ImageIO.write(grid.getImage(), "png", new File(SaveAs));
		} catch (IOException e) {
			System.err.println("Error in Writing File" + SaveAs);
		}
		
		return;
	}
	
	public void CloneSim(Parameters prms){
		
		PrintParams(outdir + File.separator + "Params4Clone.txt", prms);
		PrintUtil.printArrayToFile(outdir + File.separator + "Alive.txt", alive, 2*L, L);
		PrintUtil.printArrayToFile(outdir + File.separator + "Stress.txt", stress, L, L);
		PrintUtil.printArrayToFile(outdir + File.separator + "Residual_Stress.txt", Sr, L, L);
		PrintUtil.printArrayToFile(outdir + File.separator + "Critical_Stress.txt", Sc, L, L);
		PrintUtil.printArrayToFile(outdir + File.separator + "Stress_so_Far.txt", SsoFar, L, L);
		WriteDataHeader(outdir + File.separator + "Data_LL.txt",1);
		TakeData(outdir + File.separator + "Data_LL_1.txt",outdir + File.separator + "Data_LL_2.txt");
		
		return;
	}

		

	
	
}
package chris.foo.ofc.clusters;


import java.io.File;

import scikit.jobs.Job;
import scikit.jobs.params.Parameters;
import chris.foo.ofc.old.NfailDamage2D;
import chris.util.PrintUtil;

public class DamageClusters2D extends NfailDamage2D{
	
	// Parameters
	 public String outfileC;
	 public Boolean Cdata, Percolate;
	 
	 public ClustersV4 cluster;
	 
	 
	// Constructor
	public DamageClusters2D(Parameters params) {
		
		super(params); // I assume this is a call to SimpleDamage2D's constructor
	
		PseudoConstructorClusters(params);
		
	}
	
	
	public void PseudoConstructorClusters(Parameters params){
		
		if(params.sget("Take Cluster Data").equals("On")){
			Cdata = true;
			cluster = new ClustersV4(L,BCs);
		}
		else{
			Cdata = false;
		}

		outfileC=outdir+File.separator+"Damage_Cluster.txt";
		
		//outfileCHK=outdir+File.separator+"DamageCHK.txt";
		
		Percolate = false;
		
		return;
	}
	
	public void Avalanche() {
		
		SonFSindex=0;
		showernumber=0;
		Nshowers=0;
		
		String PLACEHOLDER = "foobar";
		
		time++;
		tkip+=dtkip;
		
		DaboutFS = FSdensity(imax);
		
		// Distribute stress from initial failure
		
		dead[0]=imax;
		search = 1;
		DistStress(PLACEHOLDER);
//		if(Percolate){
//			Job.animate();
//			System.out.println("System has percolated!");
//			return;
//		}
		
		// reset plate
		resetPlate();

		
		// search for an avalanche
		findAvalanche();
		
		while (search>0){
			
			showernumber++;
			// redistribute stress of failure(s)
			DistStress(PLACEHOLDER);
//			if(Percolate){
//				Job.animate();
//				System.out.println("System has percolated!");
//				return;
//			}
			
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
			//System.out.println("All Sites Failed!");
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

		// HowToDumpStress = evenly
		// angular bias (and then which angular bias)
		//				elipse
		//				uniform + random angle ~ cos\theta
		//				all at one cos\theta
		//				multiple cos\theta # ~ \pi R
		
		
		// if string means to distribute stress evenly . . . 
		
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
					// seems like this is where I add the site to cluster
					release=stress[dead[i]]/Nalive;
					if(Cdata) Percolate = cluster.addSite(dead[i]);
					//if(Percolate) return;
				}
				for (int j = 0; j<nbs.length; j++){
					
					if (Nbool){
						stress[nbs[j]]+=(1-alphawidth*rand.nextGaussian()-alpha)*release*alive[nbs[j]+N];
					}
					else{
						stress[nbs[j]]+=(1-alpha)*release*alive[nbs[j]+N];
					}
					
				}
			}
		}
		
		return;
	}
	
	public void WriteDataHeader(){
	
		PrintUtil.printlnToFile(outfile1,"Time","t_kip","N_avlnchs","N_dead","Rgyr","Omega","<FS_stress>","rho_FS");
		PrintUtil.printlnToFile(outfile2,"Time","t_kip","Nlives=0","Nlives=1","Nlives=2",". . .","Nlives=Nmax");
		PrintUtil.printlnToFile(outfileC,"Time","t_kip","Size of Largest Cluster");
		
		return;
	}
	
	public void WriteDataHeader(String fout, int Num){
		
		switch(Num){
		
			case 1:
				PrintUtil.printlnToFile(fout,"Time","t_kip","N_avlnchs","N_dead","Rgyr","Omega","<FS_stress>","rho_FS");
				break;
			case 2:
				PrintUtil.printlnToFile(fout,"Time","t_kip","Nlives=0","Nlives=1","Nlives=2",". . .","Nlives=Nmax");
				break;
			case 3:
				PrintUtil.printlnToFile(fout,"Time","t_kip","Size of Largest Cluster");
				break;
			default:
				System.err.println("File Not Found!");
				break;
		}
		return;
	}
	
	public void TakeData(){
		 double rgyr;
		
		int[] LS = LiveSites(); 
		int LSlength = LS.length;
		if(LSlength>0){
			rgyr=radiusGyration(LS[rand.nextInt(LSlength)]);
		}
		else{
			rgyr=0;
		}
		
		PrintUtil.printlnToFile(outfile1,time,tkip,Nshowers,NdeadS,rgyr,EFmetric(),GetAve(SonFS,SonFSindex),DaboutFS);
		PrintUtil.print2TimeAndVectorToFile(outfile2, time, tkip, NlivesLeft);
		if (Cdata) PrintUtil.printlnToFile(outfileC,time,tkip,(double)(cluster.getLargestCluster()));
		
		return;
	}
	
	public void TakeData(String fout){
		 double rgyr;
		
		int[] LS = LiveSites(); 
		int LSlength = LS.length;
		if(LSlength>0){
			rgyr=radiusGyration(LS[rand.nextInt(LSlength)]);
		}
		else{
			rgyr=0;
		}
		
		PrintUtil.printlnToFile(fout,time,tkip,Nshowers,NdeadS,rgyr,EFmetric(),GetAve(SonFS,SonFSindex),DaboutFS);
				
		return;
	}
	
}
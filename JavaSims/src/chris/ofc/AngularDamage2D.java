package chris.ofc;

import scikit.jobs.Job;
import scikit.jobs.params.Parameters;

public class AngularDamage2D extends NfailDamage2D {

	public AngularDamage2D(Parameters params) {
		super(params);
		
	}
	
	public void Avalanche() {
					
		SonFSindex=0;
		showernumber=0;
		Nshowers=0;
	
		// get density of failed sites about imax site
	
		DaboutFS = FSdensity(imax);
	
		time++;
	
		if(Nbool) alphawidth = 0;
	
	
		// redistribute stress from the failed site
	
//		int[] nbs = neighbors.get(imax);
//	
//		Nalive = 0;
//		for (int i = 0; i<nbs.length; i++){
//			Nalive+=alive[nbs[i]+N];
//		}
//		if(Nalive>0){
//			if (alive[imax]>0){
//				release=(1-alpha)*(stress[imax]-Sr[imax])/Nalive;
//			}
//			else{
//				release=(1-alpha)*stress[imax]/Nalive;
//			}
//			for (int i = 0; i<nbs.length; i++){
//				stress[nbs[i]]+=release*alive[nbs[i]+N];
//			}
//		}
	
		// reset plate conditionally
	
		if (alive[imax]>0) {
			stress[imax]=Sr[imax];
			alive[imax+N]=1;
		}
		else {
			stress[imax]=-2;
		}
	
	
		// search for the beginning of an avalanche
	
		search=0;
		for (int i = 0; i<N ; i++){
			if(stress[i]>Sc[i]){
				SonFS[SonFSindex++]=stress[i];
				dead[search++]=i;
				alive[i+N]=0;
				alive[i]--;
			}
		}
	
		Nshowers+=search;
	
		while (search>0){
			showernumber++;
			// Redistribute avalanche stress
	
			for (int i = 0; i<search; i++){
	
//				nbs = neighbors.get(dead[i]);
//				Nalive=0;
//				for (int j = 0; j<nbs.length; j++){
//					Nalive+=alive[nbs[j]+N];
//				}			
//				if(Nalive>0){			
//					if(alive[dead[i]]>0){
//						release=(1-alpha)*(stress[dead[i]]-Sr[dead[i]])/Nalive;
//					}
//					else{
//						release=(1-alpha)*stress[dead[i]]/Nalive;
//					}
//					for (int j = 0; j<nbs.length; j++){
//						stress[nbs[j]]+=release*alive[nbs[j]+N];
//					}
//				}
			}
	
			// reset plate conditionally
	
			for (int i = 0; i<search; i++){
				if (alive[dead[i]]>0) {
					stress[dead[i]]=Sr[dead[i]];
					alive[dead[i]+N]=1;
				}
				else {
					stress[dead[i]]=-2;
				}
			}
	
			// Look for subsequent avalanche
	
			search=0;
			for (int i = 0; i<N ; i++){
				if(stress[i]>Sc[i]){
					SonFS[SonFSindex++]=stress[i];
					dead[search]=i;
					alive[i+N]=0;
					alive[i]--;
					search++;
				}
			}
	
			Nshowers+=search;
			showering=true;
			Job.animate();	
		}

	}

}

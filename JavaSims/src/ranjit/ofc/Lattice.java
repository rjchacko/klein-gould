package ranjit.ofc;

public class Lattice {
    public Site site1[];
    public Bin binArray[];
    
    public int size;
    public int range;
    public int binnumber;
    public boolean boundary;
    public double alpha;
    public double maxStress;
    public double minStress;
    public double oldthreshold;
    public int iteration;
    public int maxBin;
    public int neighbors[];
    
    public int failedSites[] = null;
    
    Lattice(int size, int range, int binnumber, boolean boundary, double alpha){
		int R,N;
		int i=0;
		
		this.size=size;
		this.range=range;
		this.binnumber=binnumber;
		this.boundary=boundary;
		this.alpha=alpha;
		R=(2*this.range+1)*(2*this.range+1);
		N=size*size;
		this.neighbors= new int[R-1];
		this.failedSites=new int[size*size];
		for(i=0;i<N;i++){
		    failedSites[i]=N;
		}
    }
	
 


    void initializeLattice(){
		int i;	
		int N;
		N= size*size;
		maxStress=10000.5;
		minStress=maxStress-1;
		site1= new Site[N];
		maxBin=(int) (Math.floor(binnumber*maxStress)%binnumber)-1;
		//creating a lattice
		System.out.println("Creating lattice of size " + site1.length + ".\n");
		
		for(i=0;i<N;i++){
		    site1[i]=new Site(size,range,i,boundary);
		    site1[i].stress= (minStress+Math.random());
		}
		//Create bins 
		binArray = new Bin[binnumber];
		for(i=0;i<binnumber;i++){
		    binArray[i]=new Bin();
		}
		
		for(i=0;i<N;i++){
		    site1[i].calcBinnumber(this.binnumber);
		    site1[i].initbin=site1[i].bin;
		    binArray[site1[i].bin].add(site1[i]);
		}
	
    }

    public Site findInitiator(){
		int i=0;
		int j=0;
		Site Initiator = null;
		minStress=maxStress-1;
		
		for(i=0;i<binnumber-1;i++){
		    j=maxBin-i;
		    if(j>=0){
			if(!binArray[j].isEmpty()){
			    Initiator=binArray[j].findMax();			
			    if(Initiator.stress>minStress+0.5) {
				break;
			    }
			}	
		    }
		    else{
			if(!binArray[binnumber+j].isEmpty()){
			    Initiator=binArray[binnumber+j].findMax();			
			    if(Initiator.stress>minStress+0.5) {
				break;
			    }
			}	
		    }
		    
		}
		return Initiator;
    }
    
	public void findNPBCneighbors(Site X){

    	int i=0,j=0,x=0,y=0;
		int R;
	
		R=(2*this.range+1)*(2*this.range+1);

		for(x=X.location/this.size-this.range ; x<=X.location/this.size+this.range;x++){			
			for(y=X.location%this.size-this.range;y<=X.location%this.size+this.range;y++){				
				j=x*this.size+y;
				if(i<R-1 && j!=X.location){
					if(x<this.size && y<this.size && x>=0 && y>=0){
						this.neighbors[i]=j;
						i++;
					}
					else
					{
						this.neighbors[i]=X.location;
						i++;
					}
				}
			}
		}
	
		return;
    }
	
	public int Avalanche(Site Initiator){

		int i=0;
		int j=0;
		int k=0;
        int avsize = 0;
		int R=(2*range+1)*(2*range+1);
		int N= size*size;
		int failPointer=0;
		Site discharger=null;

		minStress=maxStress-1;
		double residual=0.25+0.5*(Math.random()-0.5);
		
		//increment avalanche size counter		
		
		findNPBCneighbors(Initiator);
		for(i=0;i<R-1;i++){		        
			j=neighbors[i];
			if(j != Initiator.location){
				
				//remove neighboring site from old bin	
				(binArray[site1[j].bin]).remove(site1[j]);
				
				//update stress of neighboring sites
				site1[j].stress=site1[j].stress+alpha*(Initiator.stress-(minStress+residual));
				
				
				//insert neighboring site in new bin
				site1[j].calcBinnumber(this.binnumber);
				site1[j].updatedAt=iteration;
				
				binArray[site1[j].bin].add(site1[j]);
				
				
				//add site to list of avalanching sites if over threshold
				if(site1[j].stress>maxStress && !site1[j].failed){
					site1[j].failed=true;
					failedSites[failPointer]=j;
					failPointer++;				
				}
			}
		}
		Initiator.stress=minStress+residual;
		Initiator.failCounter++;

		
		while(failedSites[k]!=N){
			
		    discharger=site1[failedSites[k]];
		    failedSites[k]=N;
		    k++;
		    residual=0.25+0.5*(Math.random()-0.5);
		    if(k==N) k=0;	
		    //increment avalanche size counter
		    avsize++;
		    //turn off fail flag
		    discharger.failed=false;
		    findNPBCneighbors(discharger);
		    for(i=0;i<R-1;i++){
			j=neighbors[i];			
			if(j != discharger.location){
			    
			    //remove receiving site from current bin
			    (binArray[site1[j].bin]).remove(site1[j]);
			    //update stress of neighbor
			    site1[j].stress=site1[j].stress+alpha*(discharger.stress-(minStress+residual));
			    
			    //insert receiving site in new bin
			    site1[j].calcBinnumber(this.binnumber);
			    site1[j].updatedAt=iteration;
			    (binArray[site1[j].bin]).add(site1[j]);
			    
			    //add any new failing sites to queue
			    if(site1[j].stress>maxStress && ! site1[j].failed){
				site1[j].failed=true;
				failedSites[failPointer]=j;
				failPointer++;				
			    }
			}
			if(failPointer==N) failPointer=0;
		    }
		    discharger.stress=minStress+residual;
		    discharger.failCounter++;
		    
		}
		
		return avsize;
	}
	
	public double averageStress(){
		int i=0;
		int N;
		double avgStress=0;
		N=size*size;
		
		for(i=0;i<N;i++){
			avgStress+=site1[i].stress-(maxStress-1);
		}
		avgStress=avgStress/N;
		
		return avgStress;
		
	}
	
}

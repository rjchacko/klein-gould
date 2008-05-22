package ranjit.ofc;

public class Lattice {
    public Site sites[];
    public Bin binArray[];
    
    public int size;
    public int range;
    public int binnumber;
    public boolean boundary;
    public double alpha;
    public double maxStress;
    public double minStress;
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
		R=(2*this.range+1)*(2*this.range+1)-1;
		N=size*size;
		this.neighbors= new int[R];
		this.failedSites=new int[N];
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
		sites= new Site[N];
		maxBin=(int) (Math.floor(binnumber*maxStress)%binnumber)-1;
		//creating a lattice
		System.out.println("Creating lattice of size " + sites.length + ".\n");
		
		for(i=0;i<N;i++){
		    sites[i]=new Site(size,range,i,boundary);
		    sites[i].stress= (minStress+Math.random());
		}
		//Create bins 
		binArray = new Bin[binnumber];
		for(i=0;i<binnumber;i++){
		    binArray[i]=new Bin();
		}
		
		for(i=0;i<N;i++){
		    sites[i].calcBinnumber(this.binnumber);
		    sites[i].initbin=sites[i].bin;
		    binArray[sites[i].bin].add(sites[i]);
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

    	int i=0,j=0;
		int R;
	
		R=(2*this.range+1)*(2*this.range+1);

		for(int x=X.location/this.size-this.range ; x<=X.location/this.size+this.range;x++){			
			for(int y=X.location%this.size-this.range;y<=X.location%this.size+this.range;y++){				
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
        int avsize = 0;
		int R=(2*range+1)*(2*range+1)-1;
		int N= size*size;
		int failPointer=0;
		Site discharger=null;

		minStress=maxStress-1;

		failedSites[0]=Initiator.location;
		failPointer++;
		int k=0;
		while(failedSites[k]!=N){			
		    discharger=sites[failedSites[k]];
		    failedSites[k]=N;
		    k++;
		    double residual=0.25+0.5*(Math.random()-0.5);		    	
		    //increment avalanche size counter
		    avsize++;
		    //turn off fail flag
		    discharger.failed=false;
		    findNPBCneighbors(discharger);
		    for(int i=0;i<R;i++){
				int j=neighbors[i];			
				if(j != discharger.location){				    
				    //remove receiving site from current bin
				    (binArray[sites[j].bin]).remove(sites[j]);
				    //update stress of neighbor
				    sites[j].stress=sites[j].stress+alpha*(discharger.stress-(minStress+residual));
				    
				    //insert receiving site in new bin
				    sites[j].calcBinnumber(this.binnumber);
				    sites[j].updatedAt=iteration;
				    (binArray[sites[j].bin]).add(sites[j]);
				    
				    //add any new failing sites to queue
				    if(sites[j].stress>maxStress && ! sites[j].failed){
						sites[j].failed=true;
						failedSites[failPointer]=j;
						failPointer++;				
				    }
				}
				if(failPointer==N) failPointer=0;
			}
			discharger.stress=minStress+residual;
			discharger.failCounter++;
			if(k==N) k=0;
		}
		
		return avsize;
	}
	
	public double averageStress(){
		int N;
		double avgStress=0;
		N=size*size;
		
		for(int i=0;i<N;i++){
			avgStress+=sites[i].stress-(maxStress-1);
		}
		avgStress=avgStress/N;
		
		return avgStress;
		
	}
	
}

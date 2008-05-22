package ranjit.ofc;


public class Site{
    public double stress;
    //public double threshold;
    public int location;
    public boolean failed=false;
    public int bin;
    public int initbin;
    public int updatedAt;
    public int failCounter=0;
    public int size;
    public int range;
    public boolean boundary;
    public Site nextSite;


    Site(int L, int R, int position, boolean BC){
			
		this.location=position;
		this.updatedAt=0;
		this.size=L;
		this.range=R;
		this.boundary=BC;
	
	}

    
    public void calcBinnumber(int binnumber){
		int calcBin;	
		calcBin=(int) (Math.floor(binnumber*this.stress)%binnumber);
		  
		if(calcBin>binnumber-1){
		    calcBin=binnumber-1;
		}
		if(calcBin<0){
		    calcBin=0;
		}
		this.bin=calcBin;
	
		return;
    }

}
		
	

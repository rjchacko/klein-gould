package ranjit.ofc;


public class Site implements Comparable<Object>{
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
    
    public int compareTo(Object anotherSite) throws ClassCastException {
	if (!(anotherSite instanceof Site))
	    throw new ClassCastException("A Site object expected.");
	
	double anotherSiteStress = ((Site) anotherSite).stress;  
	int comparison=0;
	if(this.stress>anotherSiteStress) comparison=1;
	if(this.stress<anotherSiteStress) comparison=-1;
	if(this.stress==anotherSiteStress) comparison=0;
	
	return comparison;
	
    }

}
		
	

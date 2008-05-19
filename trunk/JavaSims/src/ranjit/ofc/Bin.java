package ranjit.ofc;

    public class Bin {
	
    	public Site firstSite=null;	
	
	public Bin(){
	    return;

	}

	public boolean remove(Site X){
	   
	    Site prevSitePointer=null;    
	    Site SitePointer=this.firstSite;

	    //find Site X in bin
	    while(SitePointer != X){
		prevSitePointer=SitePointer;
		SitePointer= SitePointer.nextSite;
	    }
	    	    
	    if(SitePointer==X){
	    		if(SitePointer==this.firstSite){
	    			this.firstSite= this.firstSite.nextSite;
	    		}
	
	    		else if (SitePointer!=this.firstSite){
	    			prevSitePointer.nextSite=SitePointer.nextSite;
	    		}
		//set forward links to zero
		SitePointer.nextSite=null;
	
	    }
	    else{
		//fail fast
		System.out.println("Binning Failed!!!!!");
		return false;
	    }
	    X.nextSite=null;
	    return true; 

	}
	
	public boolean add(Site X){
	    
	    Site SitePointer=this.firstSite;
   
	    if(SitePointer==null){
		this.firstSite=X;
	    }
	    else if(SitePointer!=null && SitePointer.nextSite==null){
		this.firstSite.nextSite=X;
		X.nextSite=null;
	    }
	    else{	
		while(SitePointer.nextSite!=null){
		    SitePointer=SitePointer.nextSite;
		}
		SitePointer.nextSite=X;
		X.nextSite=null;
	    }

	    return true;
	}
	
	public boolean isEmpty(){
	    if(this.firstSite==null){
		return true;
	    }
	    else{
		return false;
	    }
	    
	}

	public Site findMax(){
		Site maxSite=null;
		Site checkSite=null;
		
		checkSite=this.firstSite;
		maxSite=checkSite;
		if(checkSite==null){
			System.out.println("Empty bin");
			return null;
		}
		while(checkSite.nextSite!=null){
			checkSite=checkSite.nextSite;
			if(checkSite.stress>maxSite.stress) maxSite=checkSite;
		}
		return maxSite;
		
	}

    }

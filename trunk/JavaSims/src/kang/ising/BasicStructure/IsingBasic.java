package kang.ising.BasicStructure;


import chris.util.Random;


public class IsingBasic{
	
	public IsingStructure IS;
	public Percolation IBP;
	public Random flip;
	public int Fseed;
	
	public double T;     //temperature
	public double H;     //field
	
	
	public double InitialT;
	public double InitialH;
	
	public double QuenchT;  //the temperature after the quench
	public double QuenchH;   //the temperature after the quench
	
	public double field; //field for display purpose
	public double temperature; //temperature for display purpose

	public IsingBasic(IsingStructure IS, double T, double H, int Fseed)
	{
		this.IS= IS.clone();
		this.InitialT=T;
		this.InitialH=H;
		this.T=T;
		this.H=H;
		this.Fseed= Fseed;
	    this.flip= new Random(Fseed);
	    this.IBP= new Percolation();
	}


	
	public void temperatureset(double finaltemp)
	{
		T=finaltemp;
	}
	
	public void fieldset(double finalfield)
	{
		H=finalfield;
	}
	
	public void fieldflip()
	{
	    fieldset(-H);
	}
	
	public void DoubleEnergyMetric(int spin[], int Dseed)  //calculate the energy metric for a given dilution realization by running the MCS on two systems 
	{
		
	}
	
	public void DoubleSpinMetric(int spin[], int Dseed)  //calculate the double spin metric
	{
		
	}
	

	

	

    
    public void findAFNucleation()
    {
    	
    }
    

    

	
	
	
	
}
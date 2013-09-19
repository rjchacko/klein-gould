package common.ofc;

/**
* 
*   @(#) StressLattice
*/ 

// imports
import java.util.HashSet;
import java.util.Set;
import java.util.ArrayList;


/**  
*   Basic Class for lattice of stress distributions 
*  <br>
*/
public class StressLattice{
    private ArrayList<Double> lattice;
    private double failureStress;
    int L;
    public StressLattice(int lin){
        L = lin;
        lattice = new ArrayList<Double>();
    }
    
    /**  
    *       setStress grows the lattice into the newLocation. 
    * 
    *   @param newLoc - new location to grow lattice into
    */
    public void setStress(int loc, double stress){
        lattice.set(loc,stress); 
    }

    /**  
    *       setStress2D grows the lattice into the newLocation. 
    *   Does not check if already in lattice to avoid the computation cost.
    * 
    *   @param x - x-coordinate of location to change stress for
    *   @param y - y-coordinate of location to change stress for    
    */
    public void setStress2D(int x , int y, double stress){
        lattice.set(x+y*L,stress); 
    }
    
    /**  
    *       setStress3D grows the lattice into the newLocation. 
    *   Does not check if already in lattice to avoid the computation cost.
    * 
    *   @param x - x-coordinate of location to change stress for
    *   @param y - y-coordinate of location to change stress for    
    *   @param z - z-coordinate of location to change stress for    
    */
    public void setStress3D(int x , int y, int z, double stress){
        lattice.set(x+y*L+z*L*L,stress); 
    }
    
    
    // test the class
    public static void main(String[] args) {
         
    }
}    
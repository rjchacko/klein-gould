package common.invperc;

/**
* 
*   @(#) InvCluster
*/ 

// imports
import java.util.HashSet;
import java.util.Set;


/**  
*   Basic Class for invasion percolation cluster 
*  <br>
*/
public class InvCluster{
    private Set<Integer> cluster;  // using set because dont want duplicates
    int L;
    public InvCluster(int lin){
        L = lin;
        cluster = new HashSet<Integer>();
    }
    
    /**  
    *       growCluster grows the cluster into the newLocation. 
    *   Does not check if already in cluster to avoid the computation cost.
    * 
    *   @param newLoc - new location to grow cluster into
    */
    public void growCluster(int newLoc){
        cluster.add(newLoc); 
    }

    /**  
    *       growCluster2D grows the cluster into the newLocation. 
    *   Does not check if already in cluster to avoid the computation cost.
    * 
    *   @param x - x-coordinate of new location to grow cluster into
    *   @param y - y-coordinate of new location to grow cluster into    
    */
    public void growCluster2D(int x , int y){
        cluster.add(x+y*L); 
    }
    
    /**  
    *       growCluster3D grows the cluster into the newLocation. 
    *   Does not check if already in cluster to avoid the computation cost.
    * 
    *   @param x - x-coordinate of new location to grow cluster into
    *   @param y - y-coordinate of new location to grow cluster into    
    *   @param z - z-coordinate of new location to grow cluster into    
    */
    public void growCluster3D(int x , int y, int z){
        cluster.add(x+y*L+z*L*L); 
    }
    
    /**  
    *       removeFromCluster removes the location from the cluster. 
    * 
    *   @param oldLoc - old location to remove from cluster
    */
    public void removeFromCluster(int oldLoc){
        cluster.remove((Integer)oldLoc); 
    }
    
    /**  
    *       removeFromCluster2D removes the location from the cluster. 
    * 
    *   @param x - x-coordinate of location to remove from cluster 
    *   @param y - y-coordinate of location to remove from cluster 
    */
    public void removeFromCluster2D(int x, int y){
        cluster.remove((Integer)(x+y*L)); 
    }
    
    /**  
    *       removeFromCluster3D removes the location from the cluster. 
    * 
    *   @param x - x-coordinate of location to remove from cluster 
    *   @param y - y-coordinate of location to remove from cluster 
    *   @param z - z-coordinate of location to remove from cluster 
    */
    public void removeFromCluster3D(int x, int y, int z){
        cluster.remove((Integer)(x+y*L+z*L*L)); 
    }
    
    /**  
    *       isPartOfCluster returns true if location is part of cluster.
    * 
    *   @param loc - location to check if part of cluster
    */
    public boolean isPartOfCluster(int loc){
        boolean inCluster = (cluster.contains(loc)) ?  true: false;
        return inCluster;
    }

        /**  
    *       isPartOfCluster returns true if location is part of cluster.
    * 
    *   @param x - x-coordinate of location to check if in cluster
    *   @param y - y-coordinate of location to check if in cluster
    */
    public boolean isPartOfCluster2D(int x, int y){
        boolean inCluster = (cluster.contains((x+y*L))) ?  true: false;
        return inCluster;
    }
    /**  
    *       isPartOfCluster returns true if location is part of cluster.
    * 
    *   @param x - x-coordinate of location to check if in cluster
    *   @param y - y-coordinate of location to check if in cluster
    *   @param z - z-coordinate of location to check if in cluster 
    */
    public boolean isPartOfCluster3D(int x, int y, int z){
        boolean inCluster = (cluster.contains((x+y*L+z*L*L))) ?  true: false;
        return inCluster;
    }
    
    // test the class
    public static void main(String[] args) {
        InvCluster clust = new InvCluster(10000);
        clust.growCluster(85);
        clust.growCluster(85);
        System.out.println(clust.isPartOfCluster(85));
        System.out.println(clust.isPartOfCluster(5));
        clust.removeFromCluster(85);
        System.out.println(clust.isPartOfCluster(85));
        
    }
}    
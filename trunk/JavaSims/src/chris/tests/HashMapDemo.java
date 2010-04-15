package chris.tests;

import java.util.*;

public class HashMapDemo {

@SuppressWarnings("unchecked")
public static void main(String[] args) {

    HashMap<String, Double> hm = new HashMap<String, Double>();
    hm.put("Rohit", new Double(3434.34));
    hm.put("Mohit", new Double(123.22));
    hm.put("Ashish", new Double(1200.34));
    hm.put("Khariwal", new Double(99.34));
    hm.put("Pankaj", new Double(-19.34));
    Set set = hm.entrySet();

    Iterator i = set.iterator();

    while(i.hasNext()){
      Map.Entry me = (Map.Entry)i.next();
      System.out.println(me.getKey() + " : " + me.getValue() );
    }

    //deposit into Rohit's Account
    double balance = hm.get("Rohit").doubleValue();
    hm.put("Rohit", new Double(balance + 1000));

    System.out.println("Rohit new balance : " + hm.get("Rohit"));

  }
}
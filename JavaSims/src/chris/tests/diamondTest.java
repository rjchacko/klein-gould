package chris.tests;

import chris.util.Random;
import chris.util.SortUtil;

public class diamondTest {
    public static void main(String[] args) {
    	int L = 10;
    	int LL = 10;
		int dir[] = new int[L];
		int cp[]  = new int[L];
		Random rand = new Random();

    	System.out.println("seed"+"\t"+"dx"+"\t"+"dy");
    	for (int jj = 25 ; jj < 26 ; jj++){
    		int dx   = (int)(jj > L ? Math.abs(jj-3*L)-L : jj);
    		int dy   = (int)(Math.abs(jj-2*L)-L);
    		System.out.print(jj);
    		System.out.print("\t");
    		System.out.print(dx);
    		System.out.print("\t");
    		System.out.println(dy);
       		for(int kk = 0 ; kk < Math.abs(dx) ; kk++){
    			dir[kk] = (int)(Math.signum(dx));
    		}
    		for(int kk = 0 ; kk < Math.abs(dy) ; kk++){
    			dir[kk+Math.abs(dx)] = (int)(Math.signum(dy)*LL);
    		}
       		for(int kk = 0 ; kk < L ; kk++){
       			System.out.println(dir[kk]);
       		}
   			System.out.println("___________");
    		// Shuffle the order of the steps
    		cp = SortUtil.shuffleArray(dir, rand);
     		for(int kk = 0 ; kk < L ; kk++){
       			System.out.print(dir[kk]);
        		System.out.print("\t");
       			System.out.println(cp[kk]);
       		}
    	}
    }
}


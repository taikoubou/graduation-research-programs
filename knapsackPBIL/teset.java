import jdk.javadoc.internal.doclets.toolkit.resources.doclets;

public class teset{
    public static void main(String[] args){
        int[] X = {0,1,1,1,0,0,0,1,0,0};
        int x = 0;
        double a = 1.123E-5,b=1.123E-5;

        for(int i=9;i>=0;i--){
            x = (x << 1) | X[i];
        }
        System.out.println(x);

        for(int i=0;i<1;i++){
            System.out.println("hogehgoe");
        }

        for(int i=0;i<10;i++){
            System.out.println(a + "," + b);
            a+= 1E-7;
            b+= 1E-7;
        }
    }
}

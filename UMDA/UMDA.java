import java.util.Random;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.Collections;

class UMDA{
    public static void main(final String args[]){
        final int N = 10;
        final Random rnd = new Random();
        final GeneSet set = new GeneSet(rnd);

        System.out.println("N,Parameter1,Parameter2,EmpricalParameter1,EmpricalParameter2,EmpricalParameter3,EmpricalParameter4,SelsectEmpricalParameter1,SelsectEmpricalParameter2,SelsectEmpricalParameter3,SelsectEmpricalParameter4,EliteX1,EliteX2,EliteCost");
        
        set.setSelectSet();
        set.setEmpiricalDensityfunc();
        
        set.printData(0);

        for(int i=0;i<N;i++){
            set.updateParameter();      //update parameter
            set.Sampling(rnd);      //sampling

            set.setSelectSet();
            set.setEmpiricalDensityfunc();
            set.printData(i+1);
        }
    }
}

class ProbabilityOperater{
    public double parameter1;
    public double parameter2;
    private static double alpha = 0.1;

    public ProbabilityOperater(){
        parameter1 = 0.5;
        parameter2 = 0.5;
    }

    public void UpdateParameterPBIL(final Gene Elite){
        parameter1 = parameter1 * (1.0 - alpha) + Elite.X1 * alpha;
        parameter2 = parameter2 * (1.0 - alpha) + Elite.X2 * alpha;
    }

    public void updateParameterUMDA(final ArrayList<Gene> SelectSet){
        int cnt1=0,cnt2=0;
        int N = SelectSet.size();

        for(Gene g : SelectSet){
            cnt1 += g.X1;
            cnt2 += g.X2;
        }
        parameter1 = parameter1 + alpha * ((double)cnt1/(double)N - parameter1);
        parameter2 = parameter2 + alpha * ((double)cnt2/(double)N - parameter2);
    }

    public void printParameterPBIL(){
        System.out.print("Parameter: (" + parameter1 + "," + parameter2 + ")");
    }
}

class GeneSet{
    public static int N = 10;
    public static int SN = 3;
    public ArrayList<Gene> GeneSet = new ArrayList<Gene>();
    public ArrayList<Gene> SelectSet = new ArrayList<Gene>();
    public Gene Elite = new Gene();
    public double EmpiricalDensityParameter[] = new double[4];
    public double SelectEmpricalDensityParameter[] = new double[4];
    private final ProbabilityOperater Probability = new ProbabilityOperater();

    public GeneSet(final Random rnd){     //make First Generation
        for(int i=0;i<N;i++){
            GeneSet.add(new Gene());
            GeneSet.get(i).setGene(rnd, 0.5, 0.5);
            GeneSet.get(i).setCost();
        }
        for(int i=0;i<N;i++){
            if(Elite.Cost < GeneSet.get(i).Cost)
                Elite.ObjectCopy(GeneSet.get(i));
        }
    }

    public void Sampling(final Random rnd){
        GeneSet.clear();
        SelectSet.clear();

        for(int i=0;i<N;i++){
            GeneSet.add(new Gene());
            GeneSet.get(i).setGene(rnd, Probability.parameter1, Probability.parameter2);
            GeneSet.get(i).setCost();
        }

        for(int i=0;i<N;i++){
            if(Elite.Cost < GeneSet.get(i).Cost)
                Elite.ObjectCopy(GeneSet.get(i));
        }
    }

    public void updateParameter(){
        Probability.updateParameterUMDA(SelectSet);
    }

    public void printGeneSet(){
        for(int i=0;i<N;i++){
            System.out.print(i+1 + ": (" + GeneSet.get(i).X1 + "," + GeneSet.get(i).X2 +") , Cost: " + GeneSet.get(i).Cost);
            System.out.println();
        }
        Probability.printParameterPBIL();
        System.out.println();
        System.out.println("Elite: (" + Elite.X1 + "," + Elite.X2 +") , Cost: " + Elite.Cost);
    }

    public void printParameter(){
        Probability.printParameterPBIL();
    }

    public void setEmpiricalDensityfunc(){
        final int count[] = new int[4];
        final int selectcount[] = new int[4];

        for(int i=0;i<4;i++){
            EmpiricalDensityParameter[i] = 0.0;
            SelectEmpricalDensityParameter[i] = 0.0;
        }

        for(int i=0;i<N;i++){
            count[GeneSet.get(i).getGene()]++;
        }

        for(int i=0;i<SN;i++){
            selectcount[SelectSet.get(i).getGene()]++;
        }

        for(int i=0;i<4;i++){
            EmpiricalDensityParameter[i] = (double)count[i]/(double)N;
            SelectEmpricalDensityParameter[i] =(double)selectcount[i]/(double)SN;
        }
    }
    public void printEmpricalParam(){
        System.out.print("EmpricalParameter: ");
        for(int i=0;i<4;i++)
            System.out.print(EmpiricalDensityParameter[i] + ",");
    }

    public void setSelectSet(){
        Collections.sort(GeneSet, new GeneComparator());
        Collections.reverse(GeneSet);

        SelectSet.clear();
        
        for(int i=0;i<SN;i++){
            SelectSet.add(GeneSet.get(i));
        }
    }

    public void printSelectSet(){
        System.out.println("PrintSelectSet");
        for(int i=0;i<SN;i++){
            System.out.print(i+1 + ": (" + SelectSet.get(i).X1 + "," + SelectSet.get(i).X2 +") , Cost: " + SelectSet.get(i).Cost);
            System.out.println();
        }
    }

    public void printData(int n){
        //N,parameter1,parameter2
        System.out.print(n + "," + Probability.parameter1 + "," + Probability.parameter2 + ",");

        //EmpricalParameter1...4
        for(double emparameter : EmpiricalDensityParameter)
            System.out.print(emparameter + ",");

        //SelectEmpricalParameter1...4
        for(double selEmparameter : SelectEmpricalDensityParameter)
            System.out.print(selEmparameter + ",");

        //EliteX1,EliteX2,EliteCost
        System.out.println(Elite.X1 + "," + Elite.X2 + "," + Elite.Cost);
    }
}

class Gene{
    public int X1,X2;
    public int Cost = 0;

    private int costFunc(){
        int ret = 0;

        if(this.X1 == 0 && this.X2 == 0)
            ret = 8;
        else if(this.X1 == 0 && this.X2 == 1)
            ret = 5;
        else if(this.X1 == 1 && this.X2 == 0)
            ret = 1;
        else if(this.X2 == 1 && this.X2 == 1)
            ret = 10;
        else ret = -1;

        return ret;
    }

    public void setCost(){
        this.Cost = costFunc();
    }

    public void setGene(final Random rnd, final double p1, final double p2){  //p is probability of X == 1
        final double px1 = rnd.nextDouble();
        final double px2 = rnd.nextDouble();

        if(px1 < p1)this.X1 = 1;
        else this.X1 = 0;

        if(px2 < p2)this.X2 = 1;
        else this.X2 = 0;
    }

    public int getGene(){
        return (this.X2 << 1) | (this.X1);
    }

    public void ObjectCopy(final Gene g){ //object copy method 
        this.X1 = g.X1;
        this.X2 = g.X2;
        this.Cost = g.Cost;
    }
}

class GeneComparator implements Comparator<Gene>{
    @Override
    public int compare(final Gene g1,final Gene g2){
        return g1.Cost - g2.Cost;
    }
}
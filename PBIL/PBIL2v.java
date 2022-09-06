import java.util.Random;
import java.util.ArrayList;

class PBIL2v{
    public static void main(String args[]){
        int N = 10;
        Random rnd = new Random();
        GeneSet set = new GeneSet(rnd);

        set.setEmpiricalDensityfunc();

        set.printGeneSet();
        set.printParameter();
        set.printEmpricalParam();

        System.out.println();
        /*for(int i=0;i<N;i++){
            set.updateParameter();
            System.out.print(i + 1 + ",");
            set.printParameter();
            System.out.println();
            //set.printGeneSet();
            set.Sampling(rnd);
        }
        System.out.println("終了時");
        //set.printGeneSet();
        set.printParameter();*/
        System.out.println();
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

    public void UpdateParameterPBIL(Gene Elite){
        parameter1 = parameter1 * (1.0 - alpha) + Elite.X1 * alpha;
        parameter2 = parameter2 * (1.0 - alpha) + Elite.X2 * alpha;
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
    private ProbabilityOperater Probability = new ProbabilityOperater();

    public GeneSet(Random rnd){     //make First Generation
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

    public void Sampling(Random rnd){
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
        Probability.UpdateParameterPBIL(Elite);
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
        int count[] = new int[4];
        int selectcount[] = new int[4];

        for(int i=0;i<4;i++){
            EmpiricalDensityParameter[i] = 0.0;
            //SelectEmpricalDensityParameter[i] = 0.0;
        }

        for(int i=0;i<N;i++){
            count[GeneSet.get(i).getGene()]++;
        }

        //for(int i=0;i<M;i++){
          //  selectcount[SelectSet.get(i).getGene()]++;
        //}

        for(int i=0;i<4;i++){
            EmpiricalDensityParameter[i] = (double)count[i]/(double)N;
        //    SelectEmpricalDensityParameter[i] =(double)selectcount[i]/(double)N;
        }
    }
    public void printEmpricalParam(){
        System.out.print("EmpricalParameter: ");
        for(int i=0;i<4;i++)
            System.out.print(EmpiricalDensityParameter[i] + ",");
    }
}

class Gene{
    public int X1,X2;
    public int Cost = 0;

    private int costFunc(){
        int ret = 0;

        if(this.X1 == 0 && this.X2 == 0)
            ret = 3;
        else if(this.X1 == 0 && this.X2 == 1)
            ret = 1;
        else if(this.X1 == 1 && this.X2 == 0)
            ret = 5;
        else if(this.X2 == 1 && this.X2 == 1)
            ret = 10;
        else ret = -1;

        return ret;
    }

    public void setCost(){
        this.Cost = costFunc();
    }

    public void setGene(Random rnd, double p1, double p2){  //p is probability of X == 1
        double px1 = rnd.nextDouble();
        double px2 = rnd.nextDouble();

        if(px1 < p1)this.X1 = 1;
        else this.X1 = 0;

        if(px2 < p2)this.X2 = 1;
        else this.X2 = 0;
    }

    public int getGene(){
        return (this.X2 << 1) | (this.X1);
    }

    public void ObjectCopy(Gene g){ //object copy method 
        this.X1 = g.X1;
        this.X2 = g.X2;
        this.Cost = g.Cost;
    }
}
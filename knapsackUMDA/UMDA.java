import java.util.Random;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.Collections;

//定数をまとめたクラス
class Constants{
    private Constants(){}

    public static final int LENGTH = 10;        //荷物の数
    public static final int MAX_BAG = 250;      //かばんの重量
    public static final int T = 10;            //施行する世代数

    public static final int PARAMETER_SIZE = 1024;  //rhoパラメータの数

    public static final int SELECT_N = 10;      //選択する個体数
    public static final int POP_N = 30;         //１世代の個体数
}

class UMDA{     //Main Class
    public static void main(String args[]){
        //final int N = 10;
        final Random rnd = new Random();
        final GeneSet set = new GeneSet(rnd);
        
        set.setSelectSet();
        set.setEmpiricalDensityfunc();

        //set.printSelectSet();
        //set.printData(0);
        System.out.print(set.Kullback() + ", ,");
        set.printData(0);

        for(int i=0;i<Constants.T;i++){
            set.updateParameter();      //update parameter
            set.Sampling(rnd);      //sampling

            set.setSelectSet();
            set.setEmpiricalDensityfunc();
            set.setRhoParameter();
            //set.printData(i+1);
            //set.printSelectSet();
            System.out.print(set.Kullback() + ", ,");
            set.printData(i);
        }
        System.out.println();
        System.out.println();
    }
}

//荷物のクラス
class Baggages{
    public int[] weight = {10,35,15,20,50,80,40,90,100,60};
    public int[] value =  {5 ,20,20,15,40,60,90,95,110,80};
}

class ProbabilityOperater{
    public double[] parameter;      //各遺伝子が１となる確率を表すパラメータ
    public double[] rhoparameter;
    //public double[] MarginalDistribution;
    //private static double alpha = 0.1;

    public ProbabilityOperater(){
        parameter = new double[Constants.LENGTH];
        rhoparameter = new double[Constants.PARAMETER_SIZE];
        for(int i=0;i<Constants.LENGTH;i++){
            parameter[i] = 0.5;
        }
        for(int i=0;i<Constants.PARAMETER_SIZE;i++){
            rhoparameter[i] = 1.0;
        }
        for(int i=0;i<Constants.PARAMETER_SIZE;i++){
            for(int j=0;j<Constants.LENGTH;j++){
                rhoparameter[i] *= parameter[j];
            }
        }
    }
/*
0110 -> (1-p4)*p3*p2*(1-p1)
*/
    /*public void UpdateParameterPBIL(final Gene Elite){
        parameter1 = parameter1 * (1.0 - alpha) + Elite.X1 * alpha;
        parameter2 = parameter2 * (1.0 - alpha) + Elite.X2 * alpha;
    }*/

    public void updateParameterUMDA(final ArrayList<Gene> SelectSet){       //set marginal distribution
        int[] cnt = new int[Constants.LENGTH];

        for(int j=0;j<Constants.SELECT_N;j++){
            for(int i=0;i<Constants.LENGTH;i++){
                cnt[i] += SelectSet.get(j).X[i];
            }
        }

        for(int i=0;i<Constants.LENGTH;i++){
            parameter[i] = (double)cnt[i]/(double)Constants.SELECT_N;
        }
    }

    public void setRhoParameter(){
        for(int i=0;i<Constants.PARAMETER_SIZE;i++){
            rhoparameter[i] = 1;
        }
        for(int i=0;i<Constants.PARAMETER_SIZE;i++){
            int maskbit = 0x01;
            for(int j=0;j<Constants.LENGTH;j++){
                if((i & maskbit) == maskbit)
                    rhoparameter[i] *= parameter[j];
                else
                    rhoparameter[i] *= (1-parameter[j]);
                maskbit <<= 1;
            }
        }
    }

    public double KullBack(final double[] selectRhoparameter){
        double ret = 0.0;

        for(int i=0;i<Constants.PARAMETER_SIZE;i++){
            Double SRpara = selectRhoparameter[i];
            Double Rhopara = rhoparameter[i];
            Double log = 0.0;
            //if(selectRhoparameter[i] <= 0.00001)continue;
            //double log = Math.log(selectRhoparameter[i]/rhoparameter[i]);

            if(Rhopara.compareTo(0.0) == 0)continue;

            if(SRpara.compareTo(0.00001) == -1)
                log = Math.log1p(SRpara / Rhopara);
            else
                log = Math.log(SRpara / Rhopara);

            ret += SRpara * log;
        }

        return ret;
    }

    /*public void printParameterPBIL(){
        System.out.print("Parameter: (" + parameter1 + "," + parameter2 + ")");
    }*/
}

class GeneSet{
    public ArrayList<Gene> GeneSet = new ArrayList<Gene>();
    public ArrayList<Gene> SelectSet = new ArrayList<Gene>();
    public Baggages baggages = new Baggages();
    public Gene Elite = new Gene();
    public double EmpiricalDensityParameter[] = new double[Constants.PARAMETER_SIZE];
    public double SelectEmpricalDensityParameter[] = new double[Constants.PARAMETER_SIZE];
    private final ProbabilityOperater Probability = new ProbabilityOperater();

    public GeneSet(final Random rnd){     //make First Generation
        for(int i=0;i<Constants.POP_N;i++){
            GeneSet.add(new Gene());
            GeneSet.get(i).setGene(rnd,Probability.parameter);
            GeneSet.get(i).setCost(baggages);
            GeneSet.get(i).setweight(baggages);
            if(GeneSet.get(i).weight > Constants.MAX_BAG){
                GeneSet.get(i).RandomDropGoods(rnd, baggages);
            }
        }
        for(int i=0;i<Constants.POP_N;i++){
            if(Elite.Cost < GeneSet.get(i).Cost)
                Elite.ObjectCopy(GeneSet.get(i));
        }
    }

    public void Sampling(final Random rnd){
        GeneSet.clear();
        SelectSet.clear();

        for(int i=0;i<Constants.POP_N;i++){
            GeneSet.add(new Gene());
            GeneSet.get(i).setGene(rnd, Probability.parameter);
            GeneSet.get(i).setCost(baggages);
            GeneSet.get(i).setweight(baggages);
            if(GeneSet.get(i).weight > Constants.MAX_BAG){
                GeneSet.get(i).RandomDropGoods(rnd, baggages);
            }
        }

        /*for(int i=0;i<Constants.POP_N;i++){
            if(Elite.Cost < GeneSet.get(i).Cost)
                Elite.ObjectCopy(GeneSet.get(i));
        }*/
    }

    public void updateParameter(){
        Probability.updateParameterUMDA(SelectSet);
    }

    public double Kullback(){
        return Probability.KullBack(SelectEmpricalDensityParameter);
    }

    /*public void printGeneSet(){
        for(int i=0;i<Constants.POP_N;i++){
            System.out.print(i+1 + ": (" + GeneSet.get(i).X1 + "," + GeneSet.get(i).X2 +") , Cost: " + GeneSet.get(i).Cost);
            System.out.println();
        }
        Probability.printParameterPBIL();
        System.out.println();
        System.out.println("Elite: (" + Elite.X1 + "," + Elite.X2 +") , Cost: " + Elite.Cost);
    }

    public void printParameter(){
        Probability.printParameterPBIL();
    }*/

    public void setEmpiricalDensityfunc(){      //経験分布を求める関数(SelectSetのみ)
        //final int count[] = new int[Constants.PARAMETER_SIZE];
        final int selectcount[] = new int[Constants.PARAMETER_SIZE];

        for(int i=0;i<Constants.PARAMETER_SIZE;i++){
            EmpiricalDensityParameter[i] = 0.0;
            SelectEmpricalDensityParameter[i] = 0.0;
        }

        /*for(Gene g : GeneSet){
            int x = 0;

            for(int i=Constants.LENGTH-1;i>=0;i--){
                x = (x << 1) | g.X[i];
            }

            count[x]++;
        }*/

        for(Gene g : SelectSet){
            int x = 0;

            for(int i=Constants.LENGTH-1;i>=0;i--){
                x = (x << 1) | g.X[i];
            }

            selectcount[x]++;
        }

        /*for(int scnt : selectcount)
            System.out.print(scnt + ",");
        System.out.println();*/

        for(int i=0;i<Constants.PARAMETER_SIZE;i++){
            //EmpiricalDensityParameter[i] = (double)count[i]/(double)Constants.POP_N;
            SelectEmpricalDensityParameter[i] =(double)selectcount[i]/(double)Constants.SELECT_N;
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
        
        for(int i=0;i<Constants.SELECT_N;i++){
            SelectSet.add(GeneSet.get(i));
        }

        if(Elite.Cost < GeneSet.get(0).Cost)Elite.ObjectCopy(GeneSet.get(0));
    }

    public void setRhoParameter(){
        Probability.setRhoParameter();
    }

    public void printSelectSet(){
        System.out.println("PrintSelectSet");
        for(int i=0;i<Constants.SELECT_N;i++){
            System.out.print("(");
            for(int xk : SelectSet.get(i).X)
                System.out.print(xk);
            System.out.print("),");
        }
        System.out.println();
    }

    public void printData(int n){
        //System.out.print(n + ",");
        //parameter1,parameter2
        /*for(double p : Probability.parameter){
            System.out.print(p + ",");
        }*/
        //EmpricalParameter1...4
        //for(double emparameter : EmpiricalDensityParameter)
            //System.out.print(emparameter + ",");

        //SelectEmpricalParameter1...4
        //for(double selEmparameter : SelectEmpricalDensityParameter)
            //System.out.print(selEmparameter + ",");

        //rhoparameter
        //for(double rho : Probability.rhoparameter)
            //System.out.print(rho + ",");

        //EliteX1,EliteX2,EliteCost
        for(int xi : Elite.X){
            System.out.print(xi);
        }
        System.out.println("," + Elite.Cost + "," + Elite.weight);
    }
}

class Gene{
    public int[] X;
    public int weight = 0;
    public int Cost = 0;

    public Gene(){
        X = new int[Constants.LENGTH];
    }

    private int costFunc(Baggages baggages){     //Objective function(knapsackにする)
        int ret = 0;

        for(int i=0;i<Constants.LENGTH;i++){
            if(X[i] == 1)ret += baggages.value[i];
        }

        return ret;
    }

    public void setCost(Baggages baggages){
        this.Cost = costFunc(baggages);
    }

    public void setweight(final Baggages baggages){
        this.weight = 0;
        for(int i=0;i<Constants.LENGTH;i++){
            if(X[i] == 1)weight += baggages.weight[i];
        }
    }

    public void setGene(final Random rnd, final double[] p){  //p is probability of X == 1

        double[] px = new double[Constants.LENGTH];

        for(int i=0;i<Constants.LENGTH;i++){
            px[i] = rnd.nextDouble();
            if(px[i] < p[i])this.X[i] = 1;
            else this.X[i] = 0;
        }
    }

    public void RandomDropGoods(Random rnd,final Baggages b){   //荷物がおもすぎたときに条件を満たすまでランダムに捨てる
        do{
            int k = rnd.nextInt(Constants.LENGTH);
            if(X[k] == 0)continue;
            X[k] = 0;
            weight -= b.weight[k];
            Cost -= b.value[k];
        }while(weight > Constants.MAX_BAG);

        setCost(b);
        setweight(b);
    }

    public int[] getGene(){
        return X;
    }

    public void ObjectCopy(final Gene g){ //object copy method 
        for(int i=0;i<Constants.LENGTH;i++)
            this.X[i] = g.X[i];
        
        this.Cost = g.Cost;
        this.weight = g.weight;
    }
}

class GeneComparator implements Comparator<Gene>{
    @Override
    public int compare(final Gene g1,final Gene g2){
        return g1.Cost - g2.Cost;
    }
}
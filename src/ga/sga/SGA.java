/*
 * 文件名：  AGA.java 2010-4-7
 * 描述：    SGA
 * 作者：    jianghaiyang
 * 时间：    2010-4-7
 * 版本号：  v1.0
 * 评审人:    
 * 评审时间:
 */
package ga.sga;

 /**
 * <p>Title: SGA</p>
 * <p>Description: SGA</p>
 * @author  jianghaiyang
 * @version  v1.0
 */
public abstract class SGA
{
    public static final int MAXPOP = 100;

    public static final int MAXSTRING = 64;
    
    /**
     * 前代群体
     */
    public Pp[] oldPop;
    /**
     * 当代群体
     */
    public Pp[] newPop;
    public Pp[] pl;
    
    /**
     * 群体大小[2~100]
     */
    public int popsize;
    /**
     * 染色体长度[2~64]
     */
    public int lchrom;
    public int gen;
    /**
     * 遗传世代数
     */
    public int maxgen;
    /**
     * 突变次数
     */
    public int nmutation;
    /**
     * 交叉次数
     */
    public int ncross;
    /**
     * 随机交叉点基因位
     */
    public int jcross;
    /**
     * 群体最大适应度个体序号
     */
    public int maxpp;
    /**
     * 群体最小适应度个体序号
     */
    public int minpp;
    /**
     * 随机数计数器
     */
    public int jrand;
    
    /**
     * 交叉概率[0.00~1.00]
     */
    public float pcross;
    /**
     * 突变概率[0.00~1.00]
     */
    public float pmutation;
    /**
     * 群体总适应度
     */
    public float sumfitness;
    /**
     * 群体平均适应度
     */
    public float avg;
    /**
     * 群体最大适应度
     */
    public float max;
    /**
     * 群体最小适应度
     */
    public float min;
    public float seed;
    public float [] rj = new float[MAXPOP];
    public float [] oldrand = new float[MAXPOP];

    public double coef;
    
    /**
     * 个体适应度计算.
     * @param x1 
     * @return 个体适应度
     */
    public abstract float objfunc(float x1);
    
    /**
     * 群体适应度统计.
     */
    public abstract void statistics(Pp[] pop);
    
    /**
     * 选择操作.
     * @return
     */
    public abstract int select();
    
    /**
     * 贝努利试验.
     * @param probability 概率 
     * @return [0,1]
     */
    public abstract int flip(float probability);
    
    /**
     * 交叉操作
     * @param parent1 第一条染色体
     * @param parent2 第二条染色体
     * @param k5 
     * @return 1为操作成功
     */
    public abstract int crossover(char[] parent1, char[] parent2, int k5);
    

    /**
     * 变异操作.
     * @param ch
     * @return
     */
    public abstract int mutation(char ch);
    
    /**
     * 群体更新.
     */
    public abstract void generation();
    
    /**
     * 初始化.
     */
    public abstract void initialize();

    /**
     * 数据输出.
     * @param gen
     */
    public abstract void report(int gen);
    
    /**
     * 生成初始群体.
     */
    public abstract void initpop();
    
    /**
     * 控制参数输入.
     */
    public abstract void initdata();
    
    /**
     * 初始化信息输出.
     */
    public abstract void initreport();
    
    /**
     * 译码.
     * @param pp
     * @return
     */
    public abstract float decode(char[] pp);
    
    /**
     * 随机数发生器.
     * @return 生成的随机数
     */
    public abstract float random1();
    
    /**
     * 重启动随机数发生器.
     */
    public abstract void randomize1();
}


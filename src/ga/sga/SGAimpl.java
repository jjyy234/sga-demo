/*
 * 文件名：  SGAimpl.java 2010-4-7
 * 描述：    SGAimpl
 * 作者：    jianghaiyang
 * 时间：    2010-4-7
 * 版本号：  v1.0
 * 评审人:    
 * 评审时间:
 */
package ga.sga;

import java.io.IOException;

/**
 * <p>Title: SGAimpl</p>
 * <p>Description: SGAimpl</p>
 * @author  jianghaiyang
 * @version  v1.0
 */
public class SGAimpl extends SGA
{
    /**
     * 数据输出字符串长度
     */
    private int reportLength = 0;
    
    /**
     * 交叉操作
     * @param parent1 第一条染色体
     * @param parent2 第二条染色体
     * @param k5 
     * @return 1为操作成功
     */
    public int crossover(char[] parent1, char[] parent2, int k5)
    {
        int j;
        // 是否交叉
        if (0 != flip(pcross))
        {
            jcross = (int) (Math.random() * (lchrom - 1));
            ncross = ncross + 1;
        }
        else
        {
            jcross = lchrom;
        }
        
        if (jcross != lchrom)
        {
            //随机交叉点基因位之前
            for (j = 0; j < jcross; j++)
            {
                newPop[k5].chrom[j] = (char) mutation(parent1[j]);
                newPop[k5 + 1].chrom[j] = (char) mutation(parent2[j]);
            }
            //随机交叉点基因位之后
            for (j = jcross; j < lchrom; j++)
            {
                newPop[k5].chrom[j] = (char) mutation(parent2[j]);
                newPop[k5 + 1].chrom[j] = (char) mutation(parent1[j]);
            }
        }
        else
        {
            for (j = 0; j < lchrom; j++)
            {
                newPop[k5].chrom[j] = (char) mutation(parent1[j]);
                newPop[k5 + 1].chrom[j] = (char) mutation(parent2[j]);
            }
        }
        return 1;
    }
    
    /**
     * 译码.
     * @param pp
     * @return
     */
    public float decode(char[] pp)
    {
        int j;
        float tt = 0.0F;
        float tt1 = 1.0F;
        
        for (j = lchrom - 1; j > -1; j--)
        {
            //如果基因座的值是1
            if (0 != pp[j])
            {
                tt = tt + tt1;
            }
            //向高位移一位
            tt1 = 2.0F * tt1;
        }
        tt = (float) (tt / coef);
        return tt;
    }
    
    /**
     * 贝努利试验.
     * @param probability 概率 
     * @return [0,1]
     */
    public int flip(float probability)
    {
        float ppp = (float) ((Math.random() * 20001) / 20000.0);
        return ppp <= probability ? 1 : 0;
    }
    
    /**
     * 群体更新.
     */
    public void generation()
    {
        int j;
        int mate1;
        int mate2;
        j = 0;
        do
        {
            mate1 = select();
            mate2 = select();
            //交叉操作
            // 有可能
            crossover(oldPop[mate1].chrom, oldPop[mate2].chrom, j);
            
            newPop[j].x = (float) decode(newPop[j].chrom);
            newPop[j].fitness = objfunc(newPop[j].x);
            newPop[j].parent1 = mate1;
            newPop[j].parent2 = mate2;
            newPop[j].xsite = jcross;
            
            newPop[j + 1].x = (float) decode(newPop[j + 1].chrom);
            newPop[j + 1].fitness = objfunc(newPop[j + 1].x);
            newPop[j + 1].parent1 = mate1;
            newPop[j + 1].parent2 = mate2;
            newPop[j + 1].xsite = jcross;
            
            j = j + 2;
        } while (j < popsize);
    }
    
    /**
     * 控制参数输入.
     */
    public void initdata()
    {
        byte[] b = new byte[32];
        System.out.println("----------------------------------");
        System.out.println("A Simple Generic Algorithm - SGA");
        System.out.println("     (c)X.F.SHEN JUNE 1994");
        System.out.println("      All Rights Reserved");
        System.out.println("----------------------------------");
        pause();
        System.out.println("******** SAG 数据输入与初始化 ********");
        //System.out.println("******** SAG DATA ENTRY AND INITIALIZATION ********");
        System.out.println();
        try
        {
            System.out.print("请输入群体大小[2~100] ------>");
            //System.out.print("Enter population size --------->");
            System.in.read(b);
            popsize = Integer.parseInt(new String(b).trim());
            System.out.print("请输入染色体长度[2~64] ----->");
            //System.out.print("Enter chromosome lenth -------->");
            System.in.read(b);
            lchrom = Integer.parseInt(new String(b).trim());
            System.out.print("请输入最大遗传世代数 ----------->");
            //System.out.print("Enter max.generations --------->");
            System.in.read(b);
            maxgen = Integer.parseInt(new String(b).trim());
            System.out.print("请输入交叉概率[0.00~1.00] -->");
            //System.out.print("Enter crossover probability -->");
            System.in.read(b);
            pcross = Float.parseFloat(new String(b).trim());
            System.out.print("请输入突变概率[0.00~1.00] -->");
            //System.out.print("Enter mutation probability ---->");
            System.in.read(b);
            pmutation = Float.parseFloat(new String(b).trim());
        }
        catch (IOException e)
        {
            e.printStackTrace();
        }
        //randomize();
        randomize1();
        ncross = 0;
    }
    
    /**
     * 初始化.
     */
    public void initialize()
    {
        initdata();
        System.out.println("----------------------------------");
        coef = Math.pow(2.00, lchrom) - 1.0;
        reportLength = 500 + ((42 + lchrom) * popsize);
        
        System.out.println("生成初始群体...");
        initpop();
        
        System.out.println("初始群体适应度统计...");
        statistics(oldPop);
        initreport();
    }
    
    /**
     * 生成初始群体.
     */
    public void initpop()
    {
        int j, j1;
        for (j = 0; j < popsize; j++)
        {
            for (j1 = 0; j1 < lchrom; j1++)
            {
                oldPop[j].chrom[j1] = (char) (Math.random() * 2);
            }
            oldPop[j].x = (float) decode(oldPop[j].chrom);
            oldPop[j].fitness = objfunc(oldPop[j].x);
            oldPop[j].parent1 = 0;
            oldPop[j].parent2 = 0;
            oldPop[j].xsite = 0;
        }
    }
    
    /**
     * 初始化信息输出.
     */
    public void initreport()
    {
        System.out.println("----------------------------------");
        System.out.println("A Simple Generic Algorithm - SGA");
        System.out.println("     (c)X.F.SHEN JUNE 1994");
        System.out.println("      All Rights Reserved");
        System.out.println("----------------------------------");
        //System.out.println("SGA Parameters");
        System.out.println("SGA 参数");
        System.out.println("----------------------------------");
        System.out.println();
//        System.out.println("Population size (popsize) = " + popsize);
//        System.out.println("Chromosome length (lchrom) = " + lchrom);
//        System.out.println("Crossover Probability (pcross) = " + pcross);
//        System.out.println("Mutation Probability (pmutation) = " + pmutation);
        System.out.println("群体大小   (popsize)   = " + popsize);
        System.out.println("染色体长度 (lchrom)    = " + lchrom);
        System.out.println("交叉概率   (pcross)    = " + pcross);
        System.out.println("突变概率   (pmutation) = " + pmutation);
        System.out.println("----------------------------------");
        System.out.println();
//        System.out.println("Initial Population Maximum Fitness = " + max);
//        System.out.println("Initial Population Average Fitness = " + avg);
//        System.out.println("Initial Population Minimum Fitness = " + min);
//        System.out.println("Initial Population Sum of Fitness = " + sumfitness);
        System.out.println("初始群体最大适应度 = " + max);
        System.out.println("初始群体平均适应度 = " + avg);
        System.out.println("初始群体最小适应度 = " + min);
        System.out.println("初始群体总适应度   = " + sumfitness);
        pause();
    }
    
    /**
     * 变异操作.
     * @param ch
     * @return
     */
    public int mutation(char ch)
    {
        int mutate;
        mutate = flip(pmutation);
        
        // 是否突变
        if (0 != mutate)
        {
            // 突变次数加1
            nmutation = nmutation + 1;
            if (0 != (int) ch)
            {
                ch = 0;
            }
            else
            {
                ch = 1;
            }
        }
        
        // 不突变
        if (0 != (int)ch)
        {
            return 1;
        }
        else
        {
            return 0;
        }
    }
    
    /**
     * 个体适应度计算.
     * @param x1 
     * @return 个体适应度
     */
    public float objfunc(float x1)
    {
        double y = Math.sin(2.0 * Math.PI * x1);
        return (float)(y * y);
    }
    
    /**
     * 延时.
     * @return
     */
    int pause()
    {
        int j, j1;
        int x1 = 0;
        for (j = 1; j < 25; j++)
        {
            for (j1 = 1; j1 < 2; j1++)
            {
                x1 = x1 + 1;
            }
        }
        return x1;
    }

    /**
     * 重启动随机数发生器.
     */
    public void randomize1()
    {
        int i;
        // Math.randomize();
        for (i = 0; i < lchrom; i++)
        {
            oldrand[i] = (float) ((Math.random() * 30001) / 30000.0);
        }
        jrand = 0;
    }
    
    /**
     * 随机数发生器.
     * @return 生成的随机数
     */
    public float random1()
    {
        jrand++;
        if (jrand >= lchrom)
        {
            jrand = 0;
            randomize1();
        }
        return this.oldrand[jrand];
    }
    
    /**
     * 数据输出.
     * @param gen
     */
    public void report(int gen)
    {
        int k, j;
        StringBuffer strb = new StringBuffer(reportLength);
        for (j = 0; j < 79; j++)
        {
            strb.append('*');
        }
        System.out.println();
        strb.append("\n        Population Report\n");
        strb.append("        Generation ");
        strb.append(gen);
        strb.append("\n #  parents xsite string");
        for (j = 0; j < lchrom - 5; j++)
        {
            strb.append(' ');
        }
        strb.append("x         fitness\n");
        
        for (j = 0; j < 79; j++)
        {
            strb.append('+');
        }
        strb.append('\n');
        
        for (j = 0; j < popsize; j++)
        {
            // 序号
            if (j < 10)
            {
                strb.append(' ');
            }
            strb.append(j); //序号
            
            strb.append(") (");
            if (newPop[j].parent1 < 10)
            {
                strb.append(' ');
            }
            
            // 父代序号
            strb.append(newPop[j].parent1); 
            strb.append(',');
            if (newPop[j].parent2 < 10)
            {
                strb.append(' ');
            }
            strb.append(newPop[j].parent2);
            strb.append(")  ");
            
            if (newPop[j].xsite < 10)
            {
                strb.append("  ");
            }
            else if (newPop[j].xsite < 100)
            {
                strb.append(' ');
            }
            strb.append(newPop[j].xsite) //交叉位置
                .append("  ");
            
            //个体码串
            for (k = 0; k < lchrom; k++)
            {
                strb.append((int) newPop[j].chrom[k]);
            }
            
            strb.append(' ').append(newPop[j].x);  //译码值
            
            strb.append(' ').append(newPop[j].fitness) //目标函数值
                            .append("-New.\n");
        }
        for (j = 0; j < 79; j++)
        {
            strb.append('*');
        }
        
        // strb.append("\nRESULT: GEN").append(gen)
        // .append(" AVG=").append(avg)
        // .append(" MIN=").append(min)
        // .append(" MAX=").append(max)
        // .append("\nncross=").append(ncross)
        // .append(" nmutation=").append(nmutation)
        // .append('\n');
        
        strb.append("\nRESULT: 世代数")
            .append(gen)
            .append(" 平均适应度=")
            .append(avg)
            .append(" 最小适应度=")
            .append(min)
            .append(" 最大适应度=")
            .append(max)
            .append("\n        累计交叉次数=")
            .append(ncross)
            .append(" 累计突变次数=")
            .append(nmutation)
            .append('\n');
        
        System.out.println(strb.toString());
    }
    
    /**
     * 选择操作.<br/> &nbsp;&nbsp;&nbsp;&nbsp;根据适应度随机选择一条染色体
     * 
     * @return 选择的染色体序号
     */
    public int select()
    {
        double rand1;
        // 前代群体累计适应度
        double partsum = 0.0;
        // 要选择的染色体序号
        int j = 0;
        // 生成介于[0~总适应度]的随机数
        // 相当于赌轮环上的位置
        rand1 = random1() * sumfitness;
        do
        {
            partsum = partsum + oldPop[j].fitness;
            j++;
        } while ((partsum < rand1) && (j < popsize));
        //前代群体适应度累计是否小于介于[0~总适应度]的随机数
        //并且要选择的染色体序号小于群体大小则继续下一条。
        
        return j - 1;
    }
    
    /**
     * 群体适应度统计.
     */
    public void statistics(Pp[] pop)
    {
        int j;
        sumfitness = pop[0].fitness;
        min = pop[0].fitness;
        max = pop[0].fitness;
        maxpp = 0;
        minpp = 0;
        for (j = 1; j < popsize; j++)
        {
            //统计群体总适应度
            sumfitness = sumfitness + pop[j].fitness;
            //统计群体最大适应度
            if (pop[j].fitness > max)
            {
                max = pop[j].fitness;
                maxpp = j;
            }
            //统计群体最小适应度
            if (pop[j].fitness < min)
            {
                min = pop[j].fitness;
                minpp = j;
            }
        }
        //统计群体平均适应度
        avg = sumfitness / (float) popsize;
    }
    
    /**
     * 回车继续
     */
    void gerch()
    {
        try
        {
            System.in.read();
        }
        catch (IOException e)
        {
            e.printStackTrace();
        }
    }
    
    /**
     * SGA主程序.<br/>
     *   计算[0,  2*PI]上的sin(x)平方日最大值
     * @param args
     */
    public static void main(String[] args)
    {
        SGAimpl sga = new SGAimpl();
        int gen;
        int k;
        int j;
        float oldmax;
        int oldmaxpp;
        try
        {
            sga.oldPop = new Pp[SGA.MAXPOP];
            sga.newPop = new Pp[SGA.MAXPOP];
            sga.pl = new Pp[0];
        }
        catch(Exception e)
        {
            System.out.println("memory request failed!");
            return;
        }
        
        for (k = 0; k < SGA.MAXPOP; k++)
        {
            sga.oldPop[k] = new Pp();
            sga.oldPop[k].chrom[0] = '\0';
            
            sga.newPop[k] = new Pp();
            sga.newPop[k].chrom[0] = '\0';
        }
        gen = 0;
        
        //初始化
        sga.initialize();
        // clrscr();
        sga.pl = sga.newPop;
        sga.newPop = sga.oldPop;
        sga.statistics(sga.newPop);
        sga.report(gen);
        sga.newPop = sga.pl;
        sga.gerch();
        
        do
        {
            gen = gen + 1;
            
            // 前代群体最大适应度
            oldmax = sga.max;
            
            // 前代群体最大适应度个体序号
            oldmaxpp = sga.maxpp;
            
            sga.generation();
            sga.statistics(sga.newPop);
            
            //如果当代群体最大适应度小于前代群体最大适应度
            if (sga.max < oldmax)
            {
                //当代群体最小适应度个体等于前代群体最大适应度个体
                for (j = 0; j < sga.lchrom; j++)
                {
                    sga.newPop[sga.minpp].chrom[j] =
                        sga.oldPop[oldmaxpp].chrom[j];
                }
                sga.newPop[sga.minpp].x = sga.oldPop[oldmaxpp].x;
                sga.newPop[sga.minpp].fitness = sga.oldPop[oldmaxpp].fitness;
                sga.statistics(sga.newPop);
            }
            sga.report(gen);
            // 当代群体与前代群体交换
            sga.pl = sga.oldPop;
            sga.oldPop = sga.newPop;
            sga.newPop = sga.pl;
            sga.gerch();
        } while (gen < sga.maxgen);
    }
}



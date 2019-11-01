/*
 * �ļ�����  SGAimpl.java 2010-4-7
 * ������    SGAimpl
 * ���ߣ�    jianghaiyang
 * ʱ�䣺    2010-4-7
 * �汾�ţ�  v1.0
 * ������:    
 * ����ʱ��:
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
     * ��������ַ�������
     */
    private int reportLength = 0;
    
    /**
     * �������
     * @param parent1 ��һ��Ⱦɫ��
     * @param parent2 �ڶ���Ⱦɫ��
     * @param k5 
     * @return 1Ϊ�����ɹ�
     */
    public int crossover(char[] parent1, char[] parent2, int k5)
    {
        int j;
        // �Ƿ񽻲�
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
            //�����������λ֮ǰ
            for (j = 0; j < jcross; j++)
            {
                newPop[k5].chrom[j] = (char) mutation(parent1[j]);
                newPop[k5 + 1].chrom[j] = (char) mutation(parent2[j]);
            }
            //�����������λ֮��
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
     * ����.
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
            //�����������ֵ��1
            if (0 != pp[j])
            {
                tt = tt + tt1;
            }
            //���λ��һλ
            tt1 = 2.0F * tt1;
        }
        tt = (float) (tt / coef);
        return tt;
    }
    
    /**
     * ��Ŭ������.
     * @param probability ���� 
     * @return [0,1]
     */
    public int flip(float probability)
    {
        float ppp = (float) ((Math.random() * 20001) / 20000.0);
        return ppp <= probability ? 1 : 0;
    }
    
    /**
     * Ⱥ�����.
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
            //�������
            // �п���
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
     * ���Ʋ�������.
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
        System.out.println("******** SAG �����������ʼ�� ********");
        //System.out.println("******** SAG DATA ENTRY AND INITIALIZATION ********");
        System.out.println();
        try
        {
            System.out.print("������Ⱥ���С[2~100] ------>");
            //System.out.print("Enter population size --------->");
            System.in.read(b);
            popsize = Integer.parseInt(new String(b).trim());
            System.out.print("������Ⱦɫ�峤��[2~64] ----->");
            //System.out.print("Enter chromosome lenth -------->");
            System.in.read(b);
            lchrom = Integer.parseInt(new String(b).trim());
            System.out.print("����������Ŵ������� ----------->");
            //System.out.print("Enter max.generations --------->");
            System.in.read(b);
            maxgen = Integer.parseInt(new String(b).trim());
            System.out.print("�����뽻�����[0.00~1.00] -->");
            //System.out.print("Enter crossover probability -->");
            System.in.read(b);
            pcross = Float.parseFloat(new String(b).trim());
            System.out.print("������ͻ�����[0.00~1.00] -->");
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
     * ��ʼ��.
     */
    public void initialize()
    {
        initdata();
        System.out.println("----------------------------------");
        coef = Math.pow(2.00, lchrom) - 1.0;
        reportLength = 500 + ((42 + lchrom) * popsize);
        
        System.out.println("���ɳ�ʼȺ��...");
        initpop();
        
        System.out.println("��ʼȺ����Ӧ��ͳ��...");
        statistics(oldPop);
        initreport();
    }
    
    /**
     * ���ɳ�ʼȺ��.
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
     * ��ʼ����Ϣ���.
     */
    public void initreport()
    {
        System.out.println("----------------------------------");
        System.out.println("A Simple Generic Algorithm - SGA");
        System.out.println("     (c)X.F.SHEN JUNE 1994");
        System.out.println("      All Rights Reserved");
        System.out.println("----------------------------------");
        //System.out.println("SGA Parameters");
        System.out.println("SGA ����");
        System.out.println("----------------------------------");
        System.out.println();
//        System.out.println("Population size (popsize) = " + popsize);
//        System.out.println("Chromosome length (lchrom) = " + lchrom);
//        System.out.println("Crossover Probability (pcross) = " + pcross);
//        System.out.println("Mutation Probability (pmutation) = " + pmutation);
        System.out.println("Ⱥ���С   (popsize)   = " + popsize);
        System.out.println("Ⱦɫ�峤�� (lchrom)    = " + lchrom);
        System.out.println("�������   (pcross)    = " + pcross);
        System.out.println("ͻ�����   (pmutation) = " + pmutation);
        System.out.println("----------------------------------");
        System.out.println();
//        System.out.println("Initial Population Maximum Fitness = " + max);
//        System.out.println("Initial Population Average Fitness = " + avg);
//        System.out.println("Initial Population Minimum Fitness = " + min);
//        System.out.println("Initial Population Sum of Fitness = " + sumfitness);
        System.out.println("��ʼȺ�������Ӧ�� = " + max);
        System.out.println("��ʼȺ��ƽ����Ӧ�� = " + avg);
        System.out.println("��ʼȺ����С��Ӧ�� = " + min);
        System.out.println("��ʼȺ������Ӧ��   = " + sumfitness);
        pause();
    }
    
    /**
     * �������.
     * @param ch
     * @return
     */
    public int mutation(char ch)
    {
        int mutate;
        mutate = flip(pmutation);
        
        // �Ƿ�ͻ��
        if (0 != mutate)
        {
            // ͻ�������1
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
        
        // ��ͻ��
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
     * ������Ӧ�ȼ���.
     * @param x1 
     * @return ������Ӧ��
     */
    public float objfunc(float x1)
    {
        double y = Math.sin(2.0 * Math.PI * x1);
        return (float)(y * y);
    }
    
    /**
     * ��ʱ.
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
     * �����������������.
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
     * �����������.
     * @return ���ɵ������
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
     * �������.
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
            // ���
            if (j < 10)
            {
                strb.append(' ');
            }
            strb.append(j); //���
            
            strb.append(") (");
            if (newPop[j].parent1 < 10)
            {
                strb.append(' ');
            }
            
            // �������
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
            strb.append(newPop[j].xsite) //����λ��
                .append("  ");
            
            //�����봮
            for (k = 0; k < lchrom; k++)
            {
                strb.append((int) newPop[j].chrom[k]);
            }
            
            strb.append(' ').append(newPop[j].x);  //����ֵ
            
            strb.append(' ').append(newPop[j].fitness) //Ŀ�꺯��ֵ
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
        
        strb.append("\nRESULT: ������")
            .append(gen)
            .append(" ƽ����Ӧ��=")
            .append(avg)
            .append(" ��С��Ӧ��=")
            .append(min)
            .append(" �����Ӧ��=")
            .append(max)
            .append("\n        �ۼƽ������=")
            .append(ncross)
            .append(" �ۼ�ͻ�����=")
            .append(nmutation)
            .append('\n');
        
        System.out.println(strb.toString());
    }
    
    /**
     * ѡ�����.<br/> &nbsp;&nbsp;&nbsp;&nbsp;������Ӧ�����ѡ��һ��Ⱦɫ��
     * 
     * @return ѡ���Ⱦɫ�����
     */
    public int select()
    {
        double rand1;
        // ǰ��Ⱥ���ۼ���Ӧ��
        double partsum = 0.0;
        // Ҫѡ���Ⱦɫ�����
        int j = 0;
        // ���ɽ���[0~����Ӧ��]�������
        // �൱�ڶ��ֻ��ϵ�λ��
        rand1 = random1() * sumfitness;
        do
        {
            partsum = partsum + oldPop[j].fitness;
            j++;
        } while ((partsum < rand1) && (j < popsize));
        //ǰ��Ⱥ����Ӧ���ۼ��Ƿ�С�ڽ���[0~����Ӧ��]�������
        //����Ҫѡ���Ⱦɫ�����С��Ⱥ���С�������һ����
        
        return j - 1;
    }
    
    /**
     * Ⱥ����Ӧ��ͳ��.
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
            //ͳ��Ⱥ������Ӧ��
            sumfitness = sumfitness + pop[j].fitness;
            //ͳ��Ⱥ�������Ӧ��
            if (pop[j].fitness > max)
            {
                max = pop[j].fitness;
                maxpp = j;
            }
            //ͳ��Ⱥ����С��Ӧ��
            if (pop[j].fitness < min)
            {
                min = pop[j].fitness;
                minpp = j;
            }
        }
        //ͳ��Ⱥ��ƽ����Ӧ��
        avg = sumfitness / (float) popsize;
    }
    
    /**
     * �س�����
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
     * SGA������.<br/>
     *   ����[0,  2*PI]�ϵ�sin(x)ƽ�������ֵ
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
        
        //��ʼ��
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
            
            // ǰ��Ⱥ�������Ӧ��
            oldmax = sga.max;
            
            // ǰ��Ⱥ�������Ӧ�ȸ������
            oldmaxpp = sga.maxpp;
            
            sga.generation();
            sga.statistics(sga.newPop);
            
            //�������Ⱥ�������Ӧ��С��ǰ��Ⱥ�������Ӧ��
            if (sga.max < oldmax)
            {
                //����Ⱥ����С��Ӧ�ȸ������ǰ��Ⱥ�������Ӧ�ȸ���
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
            // ����Ⱥ����ǰ��Ⱥ�彻��
            sga.pl = sga.oldPop;
            sga.oldPop = sga.newPop;
            sga.newPop = sga.pl;
            sga.gerch();
        } while (gen < sga.maxgen);
    }
}



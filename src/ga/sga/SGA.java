/*
 * �ļ�����  AGA.java 2010-4-7
 * ������    SGA
 * ���ߣ�    jianghaiyang
 * ʱ�䣺    2010-4-7
 * �汾�ţ�  v1.0
 * ������:    
 * ����ʱ��:
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
     * ǰ��Ⱥ��
     */
    public Pp[] oldPop;
    /**
     * ����Ⱥ��
     */
    public Pp[] newPop;
    public Pp[] pl;
    
    /**
     * Ⱥ���С[2~100]
     */
    public int popsize;
    /**
     * Ⱦɫ�峤��[2~64]
     */
    public int lchrom;
    public int gen;
    /**
     * �Ŵ�������
     */
    public int maxgen;
    /**
     * ͻ�����
     */
    public int nmutation;
    /**
     * �������
     */
    public int ncross;
    /**
     * �����������λ
     */
    public int jcross;
    /**
     * Ⱥ�������Ӧ�ȸ������
     */
    public int maxpp;
    /**
     * Ⱥ����С��Ӧ�ȸ������
     */
    public int minpp;
    /**
     * �����������
     */
    public int jrand;
    
    /**
     * �������[0.00~1.00]
     */
    public float pcross;
    /**
     * ͻ�����[0.00~1.00]
     */
    public float pmutation;
    /**
     * Ⱥ������Ӧ��
     */
    public float sumfitness;
    /**
     * Ⱥ��ƽ����Ӧ��
     */
    public float avg;
    /**
     * Ⱥ�������Ӧ��
     */
    public float max;
    /**
     * Ⱥ����С��Ӧ��
     */
    public float min;
    public float seed;
    public float [] rj = new float[MAXPOP];
    public float [] oldrand = new float[MAXPOP];

    public double coef;
    
    /**
     * ������Ӧ�ȼ���.
     * @param x1 
     * @return ������Ӧ��
     */
    public abstract float objfunc(float x1);
    
    /**
     * Ⱥ����Ӧ��ͳ��.
     */
    public abstract void statistics(Pp[] pop);
    
    /**
     * ѡ�����.
     * @return
     */
    public abstract int select();
    
    /**
     * ��Ŭ������.
     * @param probability ���� 
     * @return [0,1]
     */
    public abstract int flip(float probability);
    
    /**
     * �������
     * @param parent1 ��һ��Ⱦɫ��
     * @param parent2 �ڶ���Ⱦɫ��
     * @param k5 
     * @return 1Ϊ�����ɹ�
     */
    public abstract int crossover(char[] parent1, char[] parent2, int k5);
    

    /**
     * �������.
     * @param ch
     * @return
     */
    public abstract int mutation(char ch);
    
    /**
     * Ⱥ�����.
     */
    public abstract void generation();
    
    /**
     * ��ʼ��.
     */
    public abstract void initialize();

    /**
     * �������.
     * @param gen
     */
    public abstract void report(int gen);
    
    /**
     * ���ɳ�ʼȺ��.
     */
    public abstract void initpop();
    
    /**
     * ���Ʋ�������.
     */
    public abstract void initdata();
    
    /**
     * ��ʼ����Ϣ���.
     */
    public abstract void initreport();
    
    /**
     * ����.
     * @param pp
     * @return
     */
    public abstract float decode(char[] pp);
    
    /**
     * �����������.
     * @return ���ɵ������
     */
    public abstract float random1();
    
    /**
     * �����������������.
     */
    public abstract void randomize1();
}


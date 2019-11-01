/*
 * 文件名：  Pp.java 2010-4-7
 * 描述：    Pp
 * 作者：    jianghaiyang
 * 时间：    2010-4-7
 * 版本号：  v1.0
 * 评审人:    
 * 评审时间:
 */
package ga.sga;

/**
 * <p>Title: Pp</p>
 * <p>Description: Pp</p>
 * @author  jianghaiyang
 * @version  v1.0
 */
public class Pp
{
    
    /**
     * 基因序列
     */
    char [] chrom = new char[SGA.MAXSTRING];
    /**
     * 译码后的值
     */
    float x;
    /**
     * 个体适应度
     */
    float fitness;
    /**
     * 父代第一条染色体序号
     */
    int parent1;
    /**
     * 父代第二条染色体序号
     */
    int parent2;
    /**
     * 父代交叉点
     */
    int xsite;
}


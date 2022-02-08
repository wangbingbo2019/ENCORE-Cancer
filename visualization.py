# -*- coding: utf-8 -*-
"""
@author: Xianan Dong
"""
import math
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

def getData(zscore_file,fc_file,num):
    '''
    Read the LCC Z-score list and the corresponding cutoff list.
    zscore_file:    File path of LCC Z-score .
    fc_file:        File path of cutoff.
    num:            Omics category.(0, 1, 2, 3, represents Transcriptome, Methylation, Somatic mutation and CNV omics.
    '''
    all_zscore,all_fc,fcc=[],[],[]
    zscore_file = open(zscore_file,'r')

    for i in zscore_file:
        n=i.strip().split("\t")
        all_zscore.append(n)  

    fc_file =open(fc_file,'r')  
    for i in fc_file:
        n=i.strip().split("\t")  
        all_fc.append(n)  
        f = i.strip().split(',')  
        fcc.append(f)

    zz,ff,zz1,ff1,z2=[],[],[],[],[]
    
    for i in all_zscore[num]:                  
        z=i.strip().split(",")
        zz.append(z)
        
    for i in zz:
        z1=[]
        for j in i:
            z1.append(float(j))
            z2.append(float(j))
        zz1.append(z1)
        
    for i in all_fc[num]:                  
        n=i.strip().split(",")
        ff.append(n)
        
    for i in ff:
        f1=[]
        for j in i:
            f1.append(float(j))
        ff1.append(f1)
        
    return ff1,zz1

def drawCLine(fc,zscore,lineColor,nodeColor,path):
    '''
    Draw CLine.
    fc:         List of cutoff.
    zscore:     List of LCC-zscore.
    lineColor:  The color of the line.
    nodeColor:  The color of the node.
    path:       The store path.
    '''
    tt=['BLCA','BRCA','CHOL','COAD','ESCA','HNSC','KIRC','KIRP','LIHC','LUAD','LUSC','PRAD','READ','THCA','UCEC']
    size=10.5
    fig, axes = plt.subplots(nrows=3, ncols=5, figsize=(5.8,3),dpi=600,constrained_layout=True)
    for i in range(3):
        for j in range(5):

            n=int(i)*5+int(j)
            axes[i][j].plot(fc[n],zscore[n],lineColor,lineWidth=0.5)
            axes[i][j].plot(fc[n],zscore[n],nodeColor,ms=0.5)

            a=min(fc[n])
            b=math.ceil(a*10**2)/(10**2)
            c=max(fc[n])
            d=math.floor(c*10**2)/(10**2)
            axes[i][j].set_ylim([min(zscore[n])-1,max(zscore[n])+2])
            axes[i][j].set_xlim(min(fc[n]),max(fc[n]))
            x= np.linspace(b,d,2)    
            axes[i][j].set_xticks(x)
            axes[i][j].tick_params(labelsize=size)
            axes[i][j].set_title(tt[n],fontsize=size,)
    plt.savefig(path,dpi = 600)  
    
def getFittedLine(zscore):
    '''
    Get fitted line.
    zscore: LCC z-score data of cancers.
    '''
    x=[]
    for j in range(10):
        for i in range(50):
            x.append(i)
        i=0
    parameter = np.polyfit(x, zscore, 5)
    p = np.poly1d(parameter,variable='x')
    fittedLine=[]
    for i in range(50):
        fittedLine.append(p(i))
    return fittedLine

def drawUCurve(zscore,titlename,downxlim,upxlim,upylim,downylim,size,color,num,path):
    '''
    Draw UCurve.
    zscore: LCC-zscore of cancers.
    titlename: The name of title.
    downxlim: Minimum abscissa.
    upxlim: Maximum abscissa.
    upylim: Maximum ordinate value.
    downylim: Minimum ordinate value.
    size:   The size of font.
    color: The color of curve.
    num: The number of cancer.
    path: The store path.
    '''
    x=[]
    for j in range(num):
        for i in range(50):
            x.append(i)
        i=0
    
    fig, axes = plt.subplots(nrows=1, ncols=1, figsize=(3, 3),dpi=600)
    sns.regplot(x,zscore,order=5,color=color,ci=95,x_estimator=np.mean,marker='.',truncate=True,scatter=False)
    plt.scatter(x,zscore,marker='.',c=color,alpha=0.2)
    axes.set_title(titlename, fontsize=size)
    axes.set_ylabel('z-score', fontsize=size)
    axes.set_xlabel('', fontsize=size)
    axes.set_ylim([downylim, upylim])
    axes.set_xlim([downxlim, upxlim])
    axes.tick_params(labelsize=size)
    plt.subplots_adjust(top = 0.8, bottom = 0.15, right = 0.96, left = 0.3, hspace = 0, wspace = 0)
    plt.savefig(path,dpi = 600)
    
def getOrdinate(data_list):
    '''
    Get ordinate.
    data_list: LCC-zscore of cancers.
    '''
    data_list_ = []
    for i in data_list:
        for j in i:
            data_list_.append(j)
    return data_list_

if __name__ == '__main__':
    zscore_file = '.\\output\\zscore_result.txt'
    fc_file = '.\\output\\fc_result.txt'
    de_f,de_z   = getData(zscore_file,fc_file,0)
    me_f,me_z   = getData(zscore_file,fc_file,1)
    mu_f,mu_z   = getData(zscore_file,fc_file,2)
    cnv_f,cnv_z = getData(zscore_file,fc_file,3)
    
    drawCLine(de_f,de_z,'r','ro','.\\output\\Transcriptome_CLine.png')
    drawCLine(me_f,me_z,'g','go','.\\output\\Methylation_CLine.png')
    drawCLine(mu_f,mu_z,'y','yo','.\\output\\Somatic_mutation_CLine.png')
    drawCLine(cnv_f,cnv_z,'b','bo','.\\output\\CNV_CLine.png')
    
    del de_z[1],de_z[5],de_z[6],de_z[7],de_z[7]
    del me_z[8],me_z[8]
    del mu_z[2]
    del cnv_z[9]
    
    de_z_ = getOrdinate(de_z)
    me_z_ = getOrdinate(me_z)
    mu_z_ = getOrdinate(mu_z)
    cnv_z_ = getOrdinate(cnv_z)
    
    drawUCurve(de_z_,'Transcriptome',-2,52,max(de_z_)+0.5,min(de_z_)-0.5,18,'r',10,'.\\output\\Transcriptome_UCurve.png')
    drawUCurve(me_z_,'Methylation',-2,52,max(me_z_)+0.5,min(me_z_)-0.5,18,'green',13,'.\\output\\Methylation_UCurve.png')
    drawUCurve(mu_z_,'Somatic mutation',-2,52,max(mu_z_)+0.5,min(mu_z_)-0.5,18,'orange',14,'.\\output\\Somatic_mutation_UCurve.png')
    drawUCurve(cnv_z_,'CNV',-2,52,max(cnv_z_)+0.5,min(cnv_z_)-0.5,18,'b',14,'.\\output\\CNV_UCurve.png')
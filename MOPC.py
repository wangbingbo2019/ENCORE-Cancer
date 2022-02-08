# -*- coding: utf-8 -*-
"""
@author: Xianan Dong
"""
import networkx as nx
import random
import numpy as np
import os

def openPPI(filename):
    '''
    Build the network.
    filename: The file path of the human interactome.
    '''
    G2 = nx.Graph()
    a=open(filename,"r")    
    for i in a:
        n=i.strip().split("\t")
        G2.add_edge(n[0],n[1])
        G2.add_node(n[0],label=n[2])
        G2.add_node(n[1],label=n[3])
    a.close()
    
    return G2

def readData(filename):
    '''
    Read the multi-omics data of cancers.
    filename: The file path of the multi-omics data of cancer.
    '''
    genenumber=[]
    DEfc,defc,defc1=[],[],[]
    MEfc,mefc,mefc1=[],[],[]
    Mutation,mut,mut1=[],[],[]
    CNV,cnv,cnv1=[],[],[]

    leng=0
    cancer=open(filename,"r")
    next(cancer)
    for i in cancer:
        n=i.strip().split("\t")
        defc.append(float(n[9]))
        mefc.append(float(n[7]))
        mut.append(float(n[1]))  
        if (float(n[2])+float(n[3])+float(n[4])+float(n[5])+float(n[6]))!=0:
            cnvr=(float(n[2])+float(n[3])+float(n[5])+float(n[6]))/(float(n[2])+float(n[3])+float(n[4])+float(n[5])+float(n[6]))
        else:
            cnvr=0
        cnv.append(cnvr)
        for node in G2.node():
            if G2.node[node]['label']==n[0]:
                genenumber.append(node)
                break
        leng+=1
      
    yz=int(leng//4)-1
    defc1=sorted(defc,reverse=True)
    d=defc1[yz]
    for k in defc:
        if k>=float(d):
            DEfc.append(k)
        else:
            DEfc.append(0.00)
            
    mefc1=sorted(mefc,reverse=True)
    m=mefc1[yz]
    for k in mefc:
        if k>=float(m):
            MEfc.append(k)
        else:
            MEfc.append(0.00)        
        
    mut1=sorted(mut,reverse=True)
    mu=mut1[yz]
    for k in mut:
        if k>=float(mu):
            Mutation.append(k)
        else:
            Mutation.append(0.00)     
        
    cnv1=sorted(cnv,reverse=True)
    c=cnv1[yz]
    for k in cnv:
        if k>=float(c):
            CNV.append(k)
        else:
            CNV.append(0.00)    
    
    de1=list([x for x in DEfc if x!=0.00])
    me1=list([x for x in MEfc if x!=0.00])
    mu1=list([x for x in Mutation if x!=0.00])
    cn1=list([x for x in CNV if x!=0.00])
    print (len(genenumber),len(de1),len(me1),len(mu1),len(cn1))
    
    return (genenumber,DEfc,MEfc,Mutation,CNV)


def lcc_zscore(genelist,data,ran,freq,maxfc):
    '''
     Calculated the connectivity significance of genes in the network.
     genelist:  The list of genes.
     data:      The perturbed data.
     ran:       Time of random.   
     freq:      The number of calculations.
     maxfc:     The maximum value of cutoff.
    '''
    zscore,fccutoff,sizeoflcc=[],[],[]
    genesnumber=[]
    minz,alfa=[],[]
    nodata,ndata=[],[]
    geneID,ngene=[],[]
 
    for j in range(len(data)):
        for i in range(len(data[j])):
            if data[j][i] != 0:
                geneID.append(genelist[i])
                nodata.append(data[j][i])
        ndata.append(nodata)
        ngene.append(geneID)
        minz.append(min(nodata))
        alfa.append((maxfc[j]-min(nodata))/freq)
        nodata=[]
        geneID=[]
        
    for u in range(freq):
        source=[]
        cutoff=[]
        y=0
        for g in range(len(minz)):
            cutoff.append(minz[g]+alfa[g]*(u+1))
  
        for x in range(len(ngene)):
            for k in range(len(ngene[x])):
                if(ndata[y][k]>=cutoff[y]):
                    if(int(ngene[x][k])!=0):
                        source.append(ngene[x][k])
            y+=1
    
        if len(source)<=0:
            print("no nodes\n")
            continue
        genes=set()
        for j in source:
            genes.add(j)
        genesnumber.append(len(genes))
        print('#genes:',len(genes))

        g=nx.subgraph(G2,genes)
        lg=max(nx.connected_components(g),key=len)
        largest=len(lg)

        print('size of LCC:',largest)
        sizeoflcc.append(largest)

        all_genes = G2.nodes()
        l_list = []
        
        for i in range(ran):
            black_nodes = random.sample(all_genes,len(genes))
    
            g2=nx.subgraph(G2,black_nodes)
            lg2=max(nx.connected_components(g2),key=len)
    
            l_list.append(len(lg2))

        l_mean = np.mean(l_list)
        l_std  = np.std(l_list)

        if l_std == 0:
            z_score = 0
        else:
            z_score = (1.*largest - l_mean)/l_std
            
        print('z_score:',z_score)
        print('fc_cutoff:',cutoff)
        zscore.append(z_score)
        fccutoff.append(cutoff)
        
    return (zscore,fccutoff,sizeoflcc,genesnumber,ndata,minz)

def refine(zlist,fclist,size):
    '''
    Cutoff optimization.
    '''
    flag=0
    max_fcvalue=max(fclist)
    for i in range(len(zlist)):

        if size[i]>6:
            flag=0
        else:
            flag=flag+1
        if flag>=2:
            max_fcvalue=fclist[i]
            break
    return max_fcvalue

if __name__ == '__main__':
    
    G2=openPPI(r'.\\input\\newHnet-2015.txt')
    print("nodes_n:",len(G2.nodes()),"edges_n:",len(G2.edges()))

    zDE,fcDE,zME,fcME,zMu,fcMu,zCNV,fcCNV=[],[],[],[],[],[],[],[]
    lccDE,lccME,lccMu,lccCNV=[],[],[],[]
    gnDE,gnME,gnMu,gnCNV=[],[],[],[]
    
    input_folder=r'.\\input\\cancer'
    
    filenames=os.listdir(input_folder)
    for filename in filenames:
        file=filename.split('.')
        print(file[0])
        path=input_folder+'\\'+filename
        genenumber,DEfc,MEfc,Mutation,CNV=readData(path)
      
        data=[DEfc]
        max_fcvalue=[]
        for y in range(len(data)):
            maxf=[]
            maxf.append(max(data[y]))
            da=[data[y]]
            (zlist,fclist,size,gnum,ndata,minz)=lcc_zscore(genenumber,da,100,50,maxf)
            max_fcvalue.append(refine(zlist,fclist,size)[0])
            (zlist,fclist,size,gnum,ndata,minz)=lcc_zscore(genenumber,data,500,50,max_fcvalue)
        zDE.append(zlist)
        fcDE.append(fclist)
        lccDE.append(size)
        gnDE.append(gnum)
    
        data=[MEfc]
        max_fcvalue=[]
        for y in range(len(data)):
            maxf=[]
            maxf.append(max(data[y]))
            da=[data[y]]
            (zlist,fclist,size,gnum,ndata,minz)=lcc_zscore(genenumber,da,100,50,maxf)
            max_fcvalue.append(refine(zlist,fclist,size)[0]) 
            (zlist,fclist,size,gnum,ndata,minz)=lcc_zscore(genenumber,data,500,50,max_fcvalue) 
        zME.append(zlist)
        fcME.append(fclist)
        lccME.append(size)  
        gnME.append(gnum)
    
        data=[Mutation]
        max_fcvalue=[]
        for y in range(len(data)):
            maxf=[]
            maxf.append(max(data[y]))
            da=[data[y]]
            (zlist,fclist,size,gnum,ndata,minz)=lcc_zscore(genenumber,da,100,50,maxf)
            max_fcvalue.append(refine(zlist,fclist,size)[0]) 
            (zlist,fclist,size,gnum,ndata,minz)=lcc_zscore(genenumber,data,500,50,max_fcvalue) 
        zMu.append(zlist)
        fcMu.append(fclist)
        lccMu.append(size)
        gnMu.append(gnum)
    
        data=[CNV]
        max_fcvalue=[]
        for y in range(len(data)):
            maxf=[]
            maxf.append(max(data[y]))
            da=[data[y]]
            (zlist,fclist,size,gnum,ndata,minz)=lcc_zscore(genenumber,da,100,50,maxf)
            max_fcvalue.append(refine(zlist,fclist,size)[0]) 
            (zlist,fclist,size,gnum,ndata,minz)=lcc_zscore(genenumber,data,500,50,max_fcvalue) 
        zCNV.append(zlist)
        fcCNV.append(fclist)
        lccCNV.append(size)
        gnCNV.append(gnum)
    
    list_zscore_all = []
    list_zscore_all.append(zDE)
    list_zscore_all.append(zME)
    list_zscore_all.append(zMu)
    list_zscore_all.append(zCNV)
    output_file_z = open('.\\output\\zscore_result.txt','w')
    for omics in list_zscore_all:
        for cancer in omics:
            for z in cancer:
                if(cancer.index(z)==len(cancer)-1):
                    output_file_z.write(str(z))
                else:
                    output_file_z.write(str(z))
                    output_file_z.write(',')
            output_file_z.write('\t')
        output_file_z.write('\n')
    output_file_z.close()
    
    list_fc_all = []
    list_fc_all.append(fcDE)
    list_fc_all.append(fcME)
    list_fc_all.append(fcMu)
    list_fc_all.append(fcCNV)
    output_file_f = open('.\\output\\fc_result.txt','w')
    for omics in list_fc_all:
        for cancer in omics:
            for fc in cancer:
                if(cancer.index(fc)==len(cancer)-1):
                    output_file_f.write(str(fc[0]))
                else:
                    output_file_f.write(str(fc[0]))
                    output_file_f.write(',')
            output_file_f.write('\t')
        output_file_f.write('\n')
    output_file_f.close()
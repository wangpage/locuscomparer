# -*-coding:utf-8-*-
import logging
from math import log, log10

import numpy as np
import pymysql
import matplotlib.pyplot as plt


import pandas as pd
import numpy as np

from scipy.stats import zscore

#读取数据
def read_metal(in_fn,marker_col='rsid',pval_col='pval'):
    d=pd.DataFrame()
    if (isinstance(in_fn,str)):
        d=pd.read_csv(in_fn,sep='\t')
        d.rename(columns={marker_col: 'rsid'}, inplace=True)
        d.rename(columns={pval_col: 'pval'}, inplace=True)
    elif (isinstance(in_fn, pd.DataFrame)):
        d=in_fn
    else:
        logging.WARN("The argument 'in_fn' must be a string or a data.frame")
    d['logp']=-np.log10(d['pval'])
    # d['logp'] = -np.log10(zscore(d['pval'], axis=0))
    # d['logp'] = zscore(d['pval'], axis=0)
    return d
#读数据
def read_metal1(in_fn,marker_col='rsid',pval_col='pval'):
    d=pd.DataFrame()
    if (isinstance(in_fn,str)):
        d=pd.read_csv(in_fn,sep='\t')
        d.rename(columns={marker_col: 'rsid'}, inplace=True)
        d.rename(columns={pval_col: 'pval'}, inplace=True)
    elif (isinstance(in_fn, pd.DataFrame)):
        d=in_fn
    else:
        logging.WARN("The argument 'in_fn' must be a string or a data.frame")
    d['logp']=-np.log10(zscore(d['pval'],axis=0))
    return d


#查询数据库
def get_position(x, genome = "hg19"):# hg38


    connect = pymysql.connect(host='locuscompare-us-west-2a.cvocub39nnri.us-west-2.rds.amazonaws.com',  # 数据库
                              user='locuscomparer',
                              password='12345678',
                              db='locuscompare',
                              charset='utf8')  # 服务器名,账户,密码，数据库名称
    cur = connect.cursor()

    # genome = match.arg(genome)  #匹配genome中的某一个数键
    # z=""
    # for x in x['rsid']:
    #     z=z+"'"+z+"'"
    select_sqli ="select rsid, chr, pos from tkg_p3v5a_{} where rsid in ('{}')".format(genome, '\',\''.join(x['rsid']))

    cur.execute(select_sqli)
    res=cur.fetchall()
    cur.close()  # 关闭游标
    connect.close()  # 关闭数据库连接
    y=pd.merge(x,pd.DataFrame(res,columns=["rsid","chr","pos"]),how='inner')
    return y
#继续查询
def retrieve_LD(chr,snp,population):
    connect = pymysql.connect(host='locuscompare-us-west-2a.cvocub39nnri.us-west-2.rds.amazonaws.com',  # 本地数据库
                              user='locuscomparer',
                              password='12345678',
                              db='locuscompare',
                              charset='utf8')  # 服务器名,账户,密码，数据库名称
    cur = connect.cursor()


    select_sqli1 = "select SNP_A, SNP_B, R2 from tkg_p3v5a_ld_chr{}_{} where SNP_A = '{}'".format(chr,population, snp)
    cur.execute(select_sqli1)
    res1 = cur.fetchall()

    select_sqli2 = "select SNP_B as SNP_A, SNP_A as SNP_B, R2 from tkg_p3v5a_ld_chr{}_{} where SNP_B = '{}'".format(chr, population, snp)
    cur.execute(select_sqli2)
    res2 = cur.fetchall()

    res=pd.DataFrame(res1+res2,columns=[ 'SNP_A'  ,     'SNP_B'   ,    'R2'])
    cur.close()  # 关闭游标
    connect.close()  # 关闭数据库连接
    return res



#读取pval1列+pval2列值最小的行，rsid的值。
def get_lead_snp(merged, snp=None ):
    if(snp is None):
        snp=merged['rsid'][(merged['pval1']+merged['pval2']).idxmin()]
    else:
        # if (not snp in merged['rsid']):
        #     logging.WARN("%s not found in the intersection of in_fn1 and in_fn2."%(snp))
        pass
    return str(snp)

def assign_color(rsid,snp,ld):
    ld=ld[ld['SNP_A']==snp]
    ld['color']=pd.cut(ld['R2'],bins=(0,0.2,0.4,0.6,0.8,1),labels=('blue','skyblue','darkgreen','orange','red'),include_lowest=True)
    color=pd.DataFrame(rsid)
    color=pd.merge(color,ld[['SNP_B','color']],left_on='rsid',right_on='SNP_B',how='left')
    # color['color']='blue'
    color=color.fillna('blue')
    if(snp in color['rsid']):
        color['color']='purple'
    else:
        zz=pd.DataFrame([[snp,None]],columns=('rsid','color'))
        zz['color']='purple'
        color=color[["rsid","color"]].append(zz)
    res=color['color']
    res.columns='rsid'

    return res


def add_label(merged, snp):
    k=0
    merged['label']=""
    for x in merged['rsid']:

        k+=1
        if x==snp:
            merged['label'].iloc[k]=snp
    return merged

def make_scatterplot(merged, title1, title2, color, shape, size, legend = True, legend_position= ['bottomright','topright','topleft'] ,snp=""):
    # fig, ax = plt.subplots()
    plt.scatter(merged["logp1"], #.dropna(axis=0,how='any',subset="logp1").dropna(axis=0,how='any',subset="logp2")
                merged["logp2"],#.dropna(axis=0,how='any',subset="logp1").dropna(axis=0,how='any',subset="logp2")
                s=10, c=color[:-1], marker=None, cmap=None, norm=None, vmin=None, vmax=None, alpha=0.8,
                linewidths=None,  edgecolors=None, plotnonfinite=False, data=None)
    sample_prod_df = merged[merged['rsid'].isin([snp])]
    plt.scatter(sample_prod_df["logp1"],sample_prod_df["logp2"],c="red")
    plt.annotate(snp,(sample_prod_df["logp1"],sample_prod_df["logp2"]))
    #解决坐标轴、标签重叠、颜色、形状、大小
    # plt.xticks([])  # 去x坐标刻度
    # plt.yticks([])  # 去y坐标刻度
    plt.xlabel('CAD GWAS -log10(P)')
    plt.ylabel('Coronary Artery eQTL -log10(P)')
    plt.title('Coronary Artery eQTL')
    legend_box=pd.DataFrame()
    if (legend==True):
        legend_position=legend_position
        if(legend_position=='bottomright'):
            legend_box=pd.DataFrame([[0.8, 0.40],[0.8 ,0.35],[0.8 ,0.30],[0.8 ,0.25],[ 0.8 ,0.20]])
        elif(legend_position == 'topright'):
            legend_box = pd.DataFrame([[0.8, 0.80],[0.8, 0.75], [0.8, 0.70], [0.8, 0.65], [0.8, 0.60]])
        else:
            legend_box = pd.DataFrame([[0.8, 0.80],[0.8, 0.75], [0.8, 0.70], [0.8, 0.65], [0.8, 0.60]])
    # plt.colorbar()
    col=["blue", "skyblue", "darkgreen", "orange", "red"]
    # for x in range(5):
    #     plt.scatter(legend_box.loc[x,0], legend_box.loc[x,1], color=col[x])
    # plt.plot(legend_box['x'], legend_box['y'], color="black")

    plt.show()
def make_locuszoom(metal,title,chr,color,shape,size,ylab_linebreak=False,snp=""):
    # p=plt.scatter(metal['x'],metal['y'],alpha=0.8)
    plt.scatter(metal["pos"],
                metal["logp"],
                s=10, c=color[:-1], marker=None, cmap=None, norm=None, vmin=None, vmax=None, alpha=0.8,
                linewidths=None, edgecolors=None, plotnonfinite=False, data=None)
    #设置x，y轴坐标轴，颜色、形状、大小、重复坐标处理
    # zzz = range(0, 4, 1)
    # plt.xticks(zzz,('12','12.5','13','13.5'))
    plt.xlabel('chr6(Mb)')
    plt.ylabel(title+' -log10(P)')
    sample_prod_df = metal[metal['rsid'].isin([snp])]
    plt.scatter(sample_prod_df["pos"], sample_prod_df["logp"], c="red")
    plt.annotate(snp, (sample_prod_df["pos"], sample_prod_df["logp"]))
    plt.title(title)
    # return p
    plt.show()

def make_combined_plot(merged, title1, title2, ld, chr, snp = None,
                               combine = True, legend = True,
                               legend_position = ('bottomright','topright','topleft'),
                               lz_ylab_linebreak=False) :
    snp = get_lead_snp(merged, snp)
    color = assign_color(merged['rsid'], snp, ld)
    shape=pd.DataFrame(merged['rsid']==snp)
    shape.replace(False,21,inplace=True)
    shape.replace(True, 23, inplace=True)
    shape.index=merged['rsid']
    size=pd.DataFrame(merged['rsid']==snp)
    size.replace(False,2,inplace=True)
    size.replace(True, 3, inplace=True)

    size.index=merged['rsid']

    merged = add_label(merged, snp)
    make_scatterplot(merged, title1, title2, color,shape, size, legend, legend_position,snp)
    metal1 = merged[['rsid', 'logp1', 'chr', 'pos', 'label']]
    metal1.columns=['rsid', 'logp', 'chr', 'pos', 'label']
    make_locuszoom(metal1, title1, chr, color, shape, size, lz_ylab_linebreak,snp)
    metal2= merged[['rsid', 'logp2', 'chr', 'pos', 'label']]
    metal2.columns=['rsid', 'logp', 'chr', 'pos', 'label']
    make_locuszoom(metal2, title2, chr, color, shape, size, lz_ylab_linebreak,snp)
    if (combine) :
        pass
        #画图



def locuscompare(in_fn1, in_fn2, marker_col1 = "rsid", pval_col1 = "pval",
                 title1 = "eQTL",marker_col2 = "rsid", pval_col2 = "pval", title2 = "GWAS",
                 snp = None, population = "EUR", combine = True, legend = True,
                 legend_position = ['bottomright','topright','topleft'],
                 lz_ylab_linebreak = False, genome = "hg19") :# hg38
    d1=read_metal(in_fn1, marker_col1, pval_col1) #读数据tsv
    d2=read_metal(in_fn2, marker_col2, pval_col2) #读数据tsv
    merged=pd.merge(d1, d2,on="rsid",suffixes = ("1", "2"),how='left')
    merged=get_position(merged,genome)
    chr=merged['chr'].drop_duplicates(keep="first")#得到基数
    if(len(chr) !=1):
        logging.warning("There must be one and only one chromosome.")
    snp = get_lead_snp(merged, snp)
    ld = retrieve_LD(chr[0], snp, population)
    p = make_combined_plot(merged, title1, title2, ld, chr, snp, combine,
                           legend, legend_position, lz_ylab_linebreak)

    # d3 = read_metal1(in_fn1, marker_col1, pval_col1)
    # d4 = read_metal1(in_fn2, marker_col2, pval_col2)
    # merged = pd.merge(d3, d4, on="rsid", suffixes=("1", "2"), how='left')
    # merged = get_position(merged, genome)
    # chr = merged['chr'].drop_duplicates(keep="first")
    # if (len(chr) != 1):
    #     logging.warning("There must be one and only one chromosome.")
    # snp = get_lead_snp(merged, snp)
    # ld = retrieve_LD(chr[0], snp, population)
    # p = make_combined_plot(merged, title1, title2, ld, chr, snp, combine,
    #                        legend, legend_position, lz_ylab_linebreak)



    return p

locuscompare(in_fn1='gwas.tsv', in_fn2='eqtl.tsv', marker_col1 = "rsid", pval_col1 = "pval",
                 title1 = "eQTL",marker_col2 = "rsid", pval_col2 = "pval", title2 = "GWAS",
                 snp = None, population = "EUR", combine = True, legend = True,
                 legend_position = ['bottomright','topright','topleft'],
                 lz_ylab_linebreak = False, genome ="hg19")# hg38

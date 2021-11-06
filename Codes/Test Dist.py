import numpy as np
import pandas as pd
import os
#os.environ["R_HOME"] = r"C:\Program Files\R\R-4.1.1"
#os.environ['path'] += r";C:\Program Files\R\R-4.1.1\bin"
from rpy2.robjects import r
from rpy2.robjects.conversion import localconverter
from rpy2.robjects import default_converter
from matplotlib import pyplot as plt
import seaborn as sns
import rpy2.robjects.pandas2ri as rpyp


def shift_dist_per_obj(tr_obj, co_obj):
    with localconverter(default_converter + rpyp.converter):
        r("library(Seurat)")
        subset = r("subset")
        GetAssayData = r(f"function(obj) obj[['RNA']]@{assay} %>% as.data.frame()")
        UpdateObj = r("""function(obj,exp){
            obj[['RNA']]@scale.data = as.matrix(exp)
            return(obj)
        } """)
        tr_exp = GetAssayData(tr_obj)
        co_exp = GetAssayData(co_obj)

        min_tr = tr_exp.apply(np.min,axis=1).min()
        min_co = co_exp.apply(np.min,axis=1).min()

        min_value = min(min_tr,min_co)

        tr_exp +=  np.abs(min_value)
        co_exp +=  np.abs(min_value)
        return UpdateObj(tr_obj, tr_exp) , UpdateObj(co_obj, co_exp) 
        
def shift_dist(tr_exp, co_exp):
    tr_exp = tr_exp.apply(np.mean, axis=1)
    co_exp = co_exp.apply(np.mean, axis=1)
    min_value = min(tr_exp.min(), co_exp.min())
    return tr_exp + np.abs(min_value), co_exp + np.abs(min_value)


def show_cluster_dist(tr_exp, co_exp, cluster, assay, do_shift):
    if do_shift:
        tr_exp, co_exp = shift_dist(tr_exp, co_exp)
    else:
        tr_exp = tr_exp.apply(np.mean, axis=1)
        co_exp = co_exp.apply(np.mean, axis=1)

    tr_exp = tr_exp.sort_values()[:int(len(tr_exp) * 0.99)]
    co_exp = co_exp.sort_values()[:int(len(co_exp) * 0.99)]
    plt.figure()
    sns.distplot(tr_exp, kde=True, label="Tr", color="red", norm_hist=True)
    sns.distplot(co_exp, kde=True, label="Co", color="blue", norm_hist=True)
    plt.title(f"cluster{cluster}:  distribution of expression with assay {assay}")
    plt.legend(loc=1, prop={'size': 6})
    plt.savefig(f"./plotOut/dist/{cluster}-{assay.split('.')[0]}")


def ExpTables(objTr, objCo, cluster, assay="data"):
    with localconverter(default_converter + rpyp.converter):
        r("library(Seurat)")
        subset = r("subset")
        GetAssayData = r(f"function(obj) obj[['RNA']]@{assay} %>% as.data.frame()")
        Trsub = subset(objTr, idents=cluster)
        Cosub = subset(objCo, idents=cluster)

        return GetAssayData(Trsub), GetAssayData(Cosub)


def scaling_the_data(exp):
    exp["non_zero"] = exp.apply(lambda x: x.astype(bool).sum(), axis=1)
    exp = exp.loc[exp.non_zero > exp.shape[1]*0.1,:]
    exp.drop("non_zero",axis=1,inplace = True)
    exp = exp.apply(lambda x: (x - x.mean()) / x.std(), axis=1)
    fig = plt.figure()
    sns.distplot(exp.mean(axis=1))
    fig = plt.show()
    return exp



def main(assy):
 #   os.chdir("../")
    do_shift = True if assy == "scale.data" else False
    r["source"]("Codes/Ligand_Receptor_pipeline.R")
    with localconverter(default_converter + rpyp.converter):
        subset = r("function (obj,clusters) subset(obj,ident = clusters)")
        conf = pd.read_csv("config.csv")
        pathTret = list(conf.loc[conf["Var"] == "objTr", "Value"])[0]
        pathControl = list(conf.loc[conf["Var"] == "objCo", "Value"])[0]
        toVec = list(conf.loc[conf["Var"] == "toVec", "Value"])[0].split(" ")
        fromVec = list(conf.loc[conf["Var"] == "fromVec", "Value"])[0].split(" ")
        objTr = r(f"objTr = readRDS('InputObj/{pathTret}')")
        objCo = r(f"objCo = readRDS('InputObj/{pathControl}')")
        objTr = subset(objTr,toVec + fromVec)
        objTCo = subset(objCo,toVec + fromVec)

        if do_shift:
            objTr, objCo = shift_dist_per_obj(objTr, objCo)

        toVec = list(conf.loc[conf["Var"] == "toVec", "Value"])[0].split(" ")
        fromVec = list(conf.loc[conf["Var"] == "fromVec", "Value"])[0].split(" ")

        for cluster in toVec + fromVec:
            tr_exp, co_exp = ExpTables(objTr, objCo, cluster, assay)
            scaling_the_data(tr_exp)
            #show_cluster_dist(tr_exp, co_exp, cluster, assay, False)


if __name__ == "__main__":
    print(os.getcwd())
    print(os.listdir())
    #assay = input("What Assay? ")
    assay = "data"
    main(assay)

import numpy as np
import pandas as pd
import os
os.environ["R_HOME"] = r"C:\Program Files\R\R-4.1.1"
os.environ['path'] += r";C:\Program Files\R\R-4.1.1\bin"
from rpy2.robjects import r
from rpy2 import robjects as ro
from rpy2.robjects.conversion import localconverter
from rpy2.robjects import default_converter
from Codes import Tf_Graph as tfg
from Codes import neo_connection_db as neo
from Codes import utils
from Codes import pyDictToR as pyd
import rpy2.robjects.pandas2ri as rpyp
import concurrent.futures
from sklearn import metrics
import pickle
import stringdb

assay = "data"


def string_network_compere(exp):
    clusterExprsstion = exp.copy()
    clusterExprsstion["not_zero"] = clusterExprsstion.apply(lambda x: x.astype(bool).sum(), axis=1)
    clusterExprsstion = clusterExprsstion.loc[clusterExprsstion.not_zero > clusterExprsstion.shape[1] * 0.15, :]
    clusterExprsstion.drop("not_zero",axis=1,inplace = True)
    clusterExprsstion = pd.DataFrame(clusterExprsstion.mean(axis=1))
    genes = clusterExprsstion.index
    string_ids = stringdb.get_string_ids(genes)

    n = int(string_ids.shape[0]/1000)
    ppis = []
    for i in range(n):
        ppi = stringdb.get_interaction_partners (string_ids.queryItem[1000*i:1000*(i+1)])
        ppi = ppi.loc[(ppi.preferredName_A.isin(genes.str.upper())) & (ppi.preferredName_B.isin(genes.str.upper()))& (ppi.tscore > 0.7)]
        ppis.append(ppi)
    
    ppi = stringdb.get_interaction_partners (string_ids.queryItem[1000*i:])
    ppi = ppi.loc[(ppi.preferredName_A.isin(genes.str.upper())) & (ppi.preferredName_B.isin(genes.str.upper()))& (ppi.tscore > 0.7)]
    ppis.append(ppi)


    ppi = stringdb.get_interaction_partners (string_ids.queryItem)
    ppi = ppi.loc[(ppi.preferredName_A.isin(genes)) & (ppi.preferredName_B.isin(genes))& (ppi.tscore > 0.7)]
    return ppi
    

def roc_test():
    conf = pd.read_csv("config.csv")
    pathTret = list(conf.loc[conf["Var"] == "objTr", "Value"])[0]
    pathControl = list(conf.loc[conf["Var"] == "objCo", "Value"])[0]
    name_obj_Tr = list(conf.loc[conf["Var"] == "nameObjTr", "Value"])[0]
    lr = pd.read_csv("files/LigandReceptorTableMouse.tsv",sep="\t")

    with localconverter(default_converter + rpyp.converter):
        r("library(Seurat)")
        readRDS = r("readRDS")
        FindMarkers = r("function(obj,id) FindMarkers(object = obj,ident.1=id,only.pos = TRUE,slot='data',verbose = FALSE)")
        GetAssayData = r(f"function(obj) as.data.frame(obj[['RNA']]@{assay})")
        subset = r("subset")
        GetObjIdent = r("function(obj) as.character(obj@active.ident)")
        obj_tert = readRDS(f"InputObj/{pathTret}")
        obj_tert = subset(obj_tert,idents=["0","1","2","3","5","7","10","11"])
        markers = FindMarkers(obj_tert,1)
        markers = markers.loc[markers.index.isin(lr.to)]
        no_markers = lr.loc[~lr.to.isin(markers.index),"to"].sample(frac=1).drop_duplicates()
        obj_tert = subset(obj_tert,idents="1")

        pa = pd.read_csv("./files/humanProtinAction.csv")
        pi = pd.read_csv("./files/humanProtinInfo.csv")
        TfDict,pi, pa = pyd.main_py_to_R(pi, pa)
        exp = GetAssayData(obj_tert)
        recp = lr.to
        recp = pd.DataFrame(recp.loc[recp.isin(pa["Input-node Gene Symbol"])].drop_duplicates())
        recp["label"] = recp.to.apply(lambda x: 1 if x in markers.index else 0)
        recp_pos = recp.loc[recp.label == 1]
        recp_neg = recp.loc[recp.label == 0]
        recp_neg = recp_neg.sample(frac=(recp_pos.shape[0] / recp_neg.shape[0]))
        recp = pd.concat([recp_pos,recp_neg],axis=0).sample(frac=1)

        df,graph_obj = tfg.DSA_anaylsis(exp, recp.to, pa, pi, TfDict,recps_for_roc = lr.to)
        fpr, tpr, thresholds = metrics.roc_curve(df.label, df.flow, pos_label=1)
        df["label"] = df.Recp.apply(lambda x: 1 if x in markers.index else 0)
        print(metrics.auc(fpr, tpr))
        df["pvalue"] =df.Recp.apply(lambda x: graph_obj.perm_p_values[x])
        fpr, tpr, thresholds = metrics.roc_curve(df.label, df.pvalue, pos_label=0)

        return metrics.auc(fpr, tpr)
           

if __name__ == "__main__":
    aucs = []
    for i in range(1):
         aucs.append(roc_test())
    print(np.mean(aucs))
    print(aucs)

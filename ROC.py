import numpy as np
import pandas as pd
import os
from pandas.io.gbq import to_gbq
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
import matplotlib.pyplot as plt 
from sklearn.metrics import roc_curve, auc
from itertools import cycle
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
    
def draw_roc_plot(fpr,tpr,name,classes):
        
    plt.figure()
    lw = 2
    colors = cycle(["aqua", "darkorange", "cornflowerblue","navy"][:len(fpr)])
    for i, color in zip(range(len(fpr)), colors):
        plt.plot(
            fpr[i],
            tpr[i],
            color=color,
            lw=lw,
            label="ROC curve for cell-type {0} (area = {1:0.2f})".format(classes[i], auc(fpr[i],tpr[i])),
        )

    plt.plot([0, 1], [0, 1], color="navy", lw=lw, linestyle="--")
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    plt.xlabel("False Positive Rate")
    plt.ylabel("True Positive Rate")
    plt.title("ROC Curve")
    plt.legend(loc="lower right")
    plt.savefig(f"C:\\Users\\Ron\\Desktop\\LRP -manuscropt\\{name}.pdf")
    #plt.savefig(r"C:\Users\Ron\Desktop\LRP -manuscropt\name.png",dpi=300)
    plt.savefig(f"C:\\Users\\Ron\\Desktop\\LRP -manuscropt\\{name}.png",dpi=300)
    plt.show()

def roc_test(ident,fpr_flow, tpr_flow, fpr_p, tpr_p):
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
        obj_tert = subset(obj_tert,idents=["1","3","5","2"])
        markers = FindMarkers(obj_tert,ident)
        markers = markers.loc[markers.index.isin(lr.to)]
        no_markers = lr.loc[~lr.to.isin(markers.index),"to"].sample(frac=1).drop_duplicates()
        obj_tert = subset(obj_tert,idents=ident)

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
        
        df,graph_obj = tfg.DSA_anaylsis(exp, recp.to, pa, pi, TfDict)
        df["label"] = df.Recp.apply(lambda x: 1 if x in markers.index else 0)
        print(df)
        fpr, tpr, _ = roc_curve(df.label, df.DSA, pos_label=1) 
        fpr_flow.append(fpr)
        tpr_flow.append(tpr)

        #draw_roc_plot(fpr,tpr)

        print(metrics.auc(fpr, tpr))
        df["pvalue"] = -1 * np.log(df.Recp.apply(lambda x: graph_obj.perm_p_values[x] + 0.01))
        fpr, tpr, _ = metrics.roc_curve(df.label, df.pvalue, pos_label=1)
        fpr_p.append(fpr)
        tpr_p.append(tpr)
        #draw_roc_plot(fpr,tpr)

        return fpr_flow, tpr_flow, fpr_p, tpr_p

def validate_tfs(TfDict=None,imp_df=None):
    conf = pd.read_csv("config.csv")
    pathTret = list(conf.loc[conf["Var"] == "objTr", "Value"])[0]
    pathControl = list(conf.loc[conf["Var"] == "objCo", "Value"])[0]
    name_obj_Tr = list(conf.loc[conf["Var"] == "nameObjTr", "Value"])[0]
    lr = pd.read_csv("files/LigandReceptorTableMouse.tsv",sep="\t")
    dfs= []
    with localconverter(default_converter + rpyp.converter):
        r("library(Seurat)")
        readRDS = r("readRDS")
        FindMarkers = r("function(obj,id) FindMarkers(object = obj,ident.1=id, min.pct = 0.35,only.pos = FALSE,slot='data',verbose = FALSE)")
        GetAssayData = r(f"function(obj) as.data.frame(obj[['RNA']]@{assay})")
        subset = r("subset")
        GetObjIdent = r("function(obj) as.character(obj@active.ident)")
        obj = readRDS(f"InputObj/{pathTret}")
        obj = subset(obj,idents=["1","3","5","2"])
        for ident in ["1","3","5","2"]:
            markers = FindMarkers(obj,ident)
            recps = markers.loc[markers.p_val_adj <= 0.05]
            markers_pos = markers.loc[markers.avg_log2FC > 0]
            markers_neg = markers.loc[markers.avg_log2FC < 0]
            recps = lr.loc[(lr.to.isin(markers_pos.index)) & (lr["from"].isin(markers_neg.index)),"to"].drop_duplicates()

           # recps = markers.loc[(markers.avg_log2FC >=0) & (markers.index.isin(lr.to))].index
            
            obj_tert = subset(obj,idents=ident)
    
            pa = pd.read_csv("./files/humanProtinAction.csv")
            pi = pd.read_csv("./files/humanProtinInfo.csv")
            if TfDict is None:
                TfDict,pi, pa = pyd.main_py_to_R(pi, pa)
            else: 
                _,pi, pa = pyd.main_py_to_R(pi, pa)
    
            exp = GetAssayData(obj_tert)        
          
            df,graph_obj = tfg.DSA_anaylsis(exp, recps, pa, pi, TfDict,markers=markers,reduce=False,do_permutation=False)
          #  return graph_obj
            _,graph_ob_baseline = tfg.DSA_anaylsis(exp, recps, pa, pi, TfDict,markers=markers,reduce=False,do_permutation=False,wights_flag=False)
          #  imp_df = graph_obj.calculate_p_value_for_tfs()
            imp_df = graph_obj.flow_to_tfs()
            imp_df_base_line = graph_ob_baseline.flow_to_tfs()
            imp_df= imp_df.merge(imp_df_base_line, left_index=True, right_index=True)
            imp_df["flow"] = imp_df.flow_x/imp_df.flow_y
            imp_df = imp_df[["flow"]]
           # return imp_df 
            #t1 = graph_obj.capcity_network
           # imp_df = graph_obj.flow_to_all_tf()
         #   imp_df = graph_obj.calculate_p_value_for_tfs()
            tfs_scores = tfg.enrch_tfs(exp,TfDict,False)
            tf_pd = pd.DataFrame({"pvalue": tfs_scores.values()},index=tfs_scores.keys())
            tf_pd = tf_pd.merge(pd.DataFrame(imp_df),left_index=True, right_index=True,how="inner")
            tf_pd = tf_pd.loc[~tf_pd.flow.isna()]
            tf_pd["gene"] = tf_pd.index
            tf_pd = pd.merge(tf_pd,markers[["avg_log2FC","p_val_adj"]],left_on="gene",right_index=True,how="inner")
            tf_pd["label"] = tf_pd["avg_log2FC"].apply(lambda x: 1 if x > 0 else 0)          
            pos = tf_pd.loc[tf_pd.label == 1]
            neg = tf_pd.loc[tf_pd.label == 0]
            n = min(pos.shape[0] , neg.shape[0])
            pos = pos.sample(frac=(n/ pos.shape[0]))
            neg = neg.sample(frac=(n/ neg.shape[0]))
            tf_pd = pd.concat([pos,neg],axis=0).sample(frac=1)
            tf_pd.pvalue = -1*np.log(tf_pd.pvalue)
            print(ident)
            print(tf_pd.corr())
           # tf_pd = tf_pd.loc[tf_pd.flow>1]#####cutoff 0.9
            dfs.append(tf_pd)
            

        return dfs

def build_tf_dorothea():
    r("library(dorothea)")
    with localconverter(default_converter + rpyp.converter):
        tfs = r("t1 = dorothea::dorothea_mm")
        tfs = tfs.loc[(tfs.confidence != "E") & (tfs.confidence != "D")]
        tfs_dict = {tf : tfs.loc[tfs.tf == tf,"target"].to_list() for tf in tfs.tf.drop_duplicates()}
    
    return tfs_dict


if __name__ == "__main__":
    fpr_flow = []
    tpr_flow = []
    fpr_p = []
    tpr_p = []
  #  imp_df = pd.read_csv(r"./files/test.csv",index_col=0)
    df = pd.read_csv(r"./files/TfTable.csv",index_col=0)
    tf_dict = {}
    for row in df.iterrows():
        if row[1][0] not in tf_dict:
            tf_dict[row[1][0]] = [row[1][1]]
        else: 
            tf_dict[row[1][0]] += [row[1][1]]
    tf_dict = build_tf_dorothea()
    imp_df = validate_tfs(TfDict=tf_dict)
    fprs = []
    tprs = []
    for df in imp_df:
        df = df.loc[~df.flow.isna()]
        fpr, tpr, _ = roc_curve(df.label, df.flow, pos_label=1)
        fprs.append(fpr)
        tprs.append(tpr)
    draw_roc_plot(fprs,tprs,"ROC TF",["macrophages","microglia","CD4","CD8"])
    #corr = []
    #for df in imp_df:
     #       temp_df = df.copy()
      #      temp_df = temp_df.loc[~temp_df.flow.isna()]
       #     temp_df.pvalue = -1*np.log(temp_df.pvalue)
        #    corr.append(temp_df.corr().iloc[0,1])
    #draw_roc_plot(fprs,tprs)
   # for i in ["1","3","5","2"]:
    #     fpr_flow, tpr_flow, fpr_p, tpr_p = roc_test(i,fpr_flow, tpr_flow, fpr_p, tpr_p)
    
   # draw_roc_plot(fpr_flow,tpr_flow)
  #  draw_roc_plot(fpr_p,tpr_p)


import numpy as np
import pandas as pd
import os
#os.environ["R_HOME"] = r"C:\Program Files\R\R-4.1.1"
#os.environ['path'] += r";C:\Program Files\R\R-4.1.1\bin"
from rpy2.robjects import r
from rpy2 import robjects as ro
from rpy2.robjects.conversion import localconverter
from rpy2.robjects import default_converter
from Codes import pyDictToR as pyd
from Codes import Tf_Graph as tfg
from Codes import neo_connection_db as neo
from Codes import utils
import rpy2.robjects.pandas2ri as rpyp
import concurrent.futures
import warnings

pd.options.mode.chained_assignment = None  # default='warn'
warnings.filterwarnings('ignore') 

assay = "data"

def LRP(path, toName, fromName, plot_path=None, thrshold=0.1, per_claster=False):
    with localconverter(default_converter + rpyp.converter):
        print(f"{toName}_{fromName}")
        r["source"]("Codes/Ligand_Receptor_pipeline.R")
        createLRtable = r("createLRtable")
        createCircosPlots = r("createCircosPlots")
        DElegenedNcolor = r("DElegenedNcolor")
        DSAlegenedNcolor = r("DSAlegenedNcolor")
        DSA_PLOT_TSNE = r("DSA_PLOT_TSNE")
        readRDS = r("readRDS")
        subset = r("subset")
        FindMarkers = r("function(obj,id) FindMarkers(object = obj,ident.1=id,only.pos = TRUE,slot='data',verbose = FALSE)")


        ProtNamesInfo = pd.read_csv("files/humanProtinInfo.csv")
        ProtActions = pd.read_csv("files/humanProtinAction.csv")
        TfDict = pyd.main_py_to_R(ProtNamesInfo, ProtActions)[0]
      
        obj = readRDS(f"InputObj/{path}")
        obj = subset(obj, ident=[toName, fromName])
        toExpression, fromExpression = ExpTables(obj, toName, fromName)
      
        lr = createLRtable(obj, toExpression, fromExpression, fromName, toName, assay, thrshold=thrshold)
        if len(lr) != 3:
            return None

        legRet, markerallL, markerallR = lr  
        DE = DElegenedNcolor(obj, fromName, toName, legRet, markerallL, markerallR )
       
        if len(DE) == 0:
            return None

        DSA_lst = tfg.DSA_anaylsis(toExpression, legRet["Receptor"], ProtActions, ProtNamesInfo, tfs=TfDict,markers=FindMarkers(obj,toName))
            
        if len(DSA_lst) == 0:
            return None

        DSA = DSAlegenedNcolor(DSA_lst[0])
        sang_recp = [rec for rec,pvalue in DSA_lst[1].perm_p_values.items() if pvalue <= 0.05]
        DSA_Table = DSA[2]
        DSA_Table = DSA_Table.loc[DSA_Table["DSA"] < np.inf, :]
        obj = dsa_score_per_cell(obj, toExpression, DSA_Table)

        if plot_path is None:
            createCircosPlots(toExpression, legRet, DE, ro.ListVector(DSA_lst[1].Tfs_Per_Recp_in_max_flow()), DSA,
                              fromName, toName, 'r', obj,de_recptors=sang_recp)
        else:
            createCircosPlots(toExpression, legRet, DE, ro.ListVector(DSA_lst[1].Tfs_Per_Recp_in_max_flow()), DSA,
                              fromName, toName, 'r', obj,de_recptors=sang_recp, path=plot_path)


        if plot_path is not None:
            DSA_PLOT_TSNE(obj, fromName, toName, plot_path)
            del obj
        else:
            DSA_PLOT_TSNE(obj, fromName, toName)
            del obj

    r("rm(list = ls())")
    r("gc()")
    return legRet.loc[legRet.Receptor.isin(sang_recp)], DSA_Table, DSA_lst[1]


def dsa_score_per_cell(obj, toExpression, DSA_Table, sacle_factor=1):
    mdata = r("function(obj,col='DSA_SCORE') ifelse(col %in% (obj@meta.data %>% names()),1,0)")
    mdata_col = r("function(obj,col='DSA_SCORE') obj@meta.data[col]")
    add_to_meta = r("""function(obj,values,cells,col='DSA_SCORE') {
            names(values) = cells
           obj[[col]] = values
           return (obj)
            }""")
    with localconverter(default_converter + rpyp.converter):
        toExpression = toExpression.loc[list(filter(lambda x: x in DSA_Table.Recp, toExpression.index))]
        for row in range(toExpression.shape[0]):
            toExpression.iloc[row, :] = (np.log1p(toExpression.iloc[row, :]) * float(
                DSA_Table.loc[DSA_Table.Recp == toExpression.index[row], :].DSA))

        if not mdata(obj)[0]:
            dsa_score = pd.Series(None, toExpression.columns)
        else:
            dsa_score = pd.Series(mdata_col(obj), toExpression.columns)

        dsa_score = pd.DataFrame(dsa_score)
        dsa_score.columns = ["DSA_SCORE"]
        dsa_per_cell = toExpression.apply(sum)
        dsa_per_cell = dsa_per_cell * sacle_factor
        dsa_per_cell = pd.DataFrame(dsa_per_cell)
        dsa_per_cell.columns = ["dsa_per_cell"]
        DSA_temp = pd.merge(dsa_score, dsa_per_cell, how="left", left_index=True, right_index=True)
        DSA_temp["DSA_SCORE"] = DSA_temp.apply(
            lambda x: x["dsa_per_cell"] if np.isnan(x["DSA_SCORE"]) else x["DSA_SCORE"], axis=1)
        return add_to_meta(obj, DSA_temp["DSA_SCORE"], DSA_temp["DSA_SCORE"].index)


def _helper_LRP(x):
    if len(x) == 3:
        return LRP(x[0], x[1], x[2])
    elif len(x) == 4:
        return LRP(x[0], x[1], x[2], plot_path=x[3])


def ExpTables(obj, toName, fromName, assay=assay):
    with localconverter(default_converter + rpyp.converter):
        r("library(Seurat)")
        subset = r("subset")
        GetAssayData = r(f"function(obj) obj[['RNA']]@{assay} %>% as.data.frame()")
        tosub = subset(obj, idents=toName)
        fromsub = subset(obj, idents=fromName)
        toExpression = GetAssayData(tosub)
        fromExpression = GetAssayData(fromsub)

        return toExpression, fromExpression


def run_pipline(args, objName, max_workers=3):
    r["source"]("Codes/Ligand_Receptor_pipeline.R")
    legRetList = {}
    DSA_Tables = {}
    DSA_Graphs = tfg.Obj_dict()
    DSA_mean = {} 

    with concurrent.futures.ProcessPoolExecutor(max_workers=max_workers) as executor:
        result = list(executor.map(_helper_LRP, args))

    with localconverter(default_converter + rpyp.converter):
        obj = r(f"readRDS('InputObj/{args[0][0]}')")
        for index, arg in enumerate(args):
            if arg[1] not in legRetList:
                legRetList[arg[1]] = {}
                DSA_Tables[arg[1]] = {}
                DSA_Graphs.add_obj(arg[1], tfg.Obj_dict())
                DSA_mean[arg[1]] = {}

            if result[index] is not None:
                legRetList[arg[1]][arg[2]] = result[index][0]
                DSA_Tables[arg[1]][arg[2]] = result[index][1]
                DSA_Graphs.return_GD_from_key(arg[1]).add_obj(arg[2], result[index][2])
                toExpression, fromExpression = ExpTables(obj, arg[1], arg[2])
                DSA_mean[arg[1]][arg[2]] = [
                    np.mean(utils.dsa_with_lig(toExpression, fromExpression, result[index][1], result[index][0]))]

        DSA_Mean = {key: neo.dict_to_table(value, "Cluster", "DSA_Mean") for key, value in DSA_mean.items()}
        utils.save_obj(legRetList, f"./outputObj/legRetLists_{objName}")
        utils.save_obj(DSA_Tables, f"./outputObj/DSA_Tables_{objName}")
        utils.save_obj(DSA_Graphs, f"./outputObj/DSA_Graphs_{objName}")
        utils.save_obj(DSA_Mean, f"./outputObj/DSA_mean_{objName}")

    r("rm(list = ls())")
    r("gc()")
    return legRetList, DSA_Tables, DSA_Graphs, DSA_Mean


def main():
    try:
        conf = pd.read_csv(r"./config.csv")
    except:
        raise Exception("no input CSV found")

    path = list(conf.loc[conf["Var"] == "obj_name", "Value"])[0]
    plotpath = list(conf.loc[conf["Var"] == "Plotpath", "Value"])[0]
    toVec = list(conf.loc[conf["Var"] == "toVec", "Value"])[0].split(" ")
    fromVec = list(conf.loc[conf["Var"] == "fromVec", "Value"])[0].split(" ")
    args = [(path, toName, fromName, plotpath) for toName in toVec for fromName in fromVec if
            fromName != toName]

    obj_name = list(conf.loc[conf["Var"] == "Project_name", "Value"])[0]
    lst_obj = run_pipline(args, obj_name)
    
    return lst_obj


def format_dsa_file(x):
    annot = pd.read_csv("Annot.csv")
    try:
        type_cluster = list(annot.loc[annot['Cluster'] == int(x), 'Type'])[0]
        subtype_cluster = list(annot.loc[annot['Cluster'] == int(x), 'subType'])[0]
    except:
        raise Exception(f"cluster {x} not found")

    if subtype_cluster is np.nan:
        return f"{x}/{type_cluster}"
    return f"{x}/{type_cluster}/{subtype_cluster}"

if __name__ == "__main__":
    print(os.getcwd())
    main()

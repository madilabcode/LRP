import numpy as np
import pandas as pd
# from bounded_pool_executor import BoundedProcessPoolExecutor
# import anndata2ri as ri
import os

os.environ["R_HOME"] = r"C:\Program Files\R\R-4.1.1"
os.environ['path'] += r";C:\Program Files\R\R-4.1.1\bin"
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
import pickle


def LRP(path, toName, fromName, plot_path=None, thrshold=0.1, per_claster=False):
    print(f"{toName}_{fromName}")
    r["source"]("Codes/Ligand_Receptor_pipeline.R")
    ProtNamesInfo = pd.read_csv("files/humanProtinInfo.csv")
    ProtActions = pd.read_csv("files/humanProtinAction.csv")
    subset = r("subset")

    with localconverter(default_converter + rpyp.converter):
        lst = pyd.main_py_to_R(ProtNamesInfo, ProtActions)
        TfDict = lst[0]
        obj = r(f"readRDS('{path}')")
        createLRtable = r("createLRtable")
        DElegenedNcolor = r("DElegenedNcolor")

        obj = subset(obj, ident=[toName, fromName])
        r("gc()")
        toExpression, fromExpression = ExpTables(obj, toName, fromName)
        try:
            LR = createLRtable(obj, toExpression, fromExpression, fromName, toName, 'counts', thrshold=thrshold)
        except Exception as e:
            LR = None
            print(e)

        if len(LR) != 3:
            del obj
            return None

        DE = DElegenedNcolor(obj, fromName, toName, LR[0], LR[1], LR[2])

        if len(DE) == 0:
            del obj
            return None

        legRet = LR[0]
        DSA_lst = tfg.DSA_anaylsis(toExpression, legRet["Receptor"], ProtActions, ProtNamesInfo, tfs=TfDict)

        if len(DSA_lst) == 0:
            del obj
            return None

        DSAlegenedNcolor = r("DSAlegenedNcolor")
        DSA = DSAlegenedNcolor(DSA_lst[0])

        DSA_Table = DSA[2]
        DSA_Table = DSA_Table.loc[DSA_Table["DSA"] < np.inf, :]

        createCircosPlots = r("createCircosPlots")
        if plot_path is None:
            createCircosPlots(toExpression, legRet, DE, ro.ListVector(DSA_lst[1].Tfs_Per_Recp_in_max_flow()), DSA,
                              fromName, toName, 'r', obj)
        else:
            createCircosPlots(toExpression, legRet, DE, ro.ListVector(DSA_lst[1].Tfs_Per_Recp_in_max_flow()), DSA,
                              fromName, toName, 'r', obj, path=plot_path)

        DSA_PLOT_TSNE = r("DSA_PLOT_TSNE")

        if per_claster:
            dsa_score_per_cell_all_cluster = r("dsa_score_per_cell_all_cluster")
            obj = dsa_score_per_cell_all_cluster(obj, toExpression, DSA_Table)

        else:
            obj = dsa_score_per_cell(obj, toExpression, DSA_Table)

        if plot_path is not None:
            DSA_PLOT_TSNE(obj, fromName, toName, plot_path)
            del obj
        else:
            DSA_PLOT_TSNE(obj, fromName, toName)
            del obj

    r("rm(list = ls())")
    r("gc()")
    return legRet, DSA_Table, DSA_lst[1]


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


def _helper_Classfier(args):
    return build_or_use_classfier(*args)


def ExpTables(obj, toName, fromName, assay="data"):
    with localconverter(default_converter + rpyp.converter):
        r("library(Seurat)")
        subset = r("subset")
        if assay == "data":
            GetAssayData = r("function(obj) obj[['RNA']]@data %>% as.data.frame()")
        else:
            GetAssayData = r(f"function(obj) GetAssayData(obj,{assay}) %>% as.data.frame()")


            GetAssayData = r("GetAssayData")
        rDataFrame = r("as.data.frame")

        tosub = subset(obj, idents=toName)
        fromsub = subset(obj, idents=fromName)

        toExpression = GetAssayData(tosub)
        fromExpression = GetAssayData(fromsub)

        return toExpression, fromExpression


def run_pipline(args, objName, max_workers=20):
    r["source"]("Codes/Ligand_Receptor_pipeline.R")
    legRetList = {}
    DSA_Tables = {}
    DSA_Graphs = tfg.Obj_dict()
    DSA_mean = {}

    # with BoundedProcessPoolExecutor(max_workers=max_workers) as executor:
    with concurrent.futures.ProcessPoolExecutor(max_workers=max_workers) as executor:
        result = list(executor.map(_helper_LRP, args))

    # result = list(filter(lambda x: x is not None, result))

    with localconverter(default_converter + rpyp.converter):
        obj = r(f"readRDS('{args[0][0]}')")
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
        tfg.save_obj(legRetList, f"./outputObj/legRetLists_{objName}")
        tfg.save_obj(DSA_Tables, f"./outputObj/DSA_Tables_{objName}")
        tfg.save_obj(DSA_Graphs, f"./outputObj/DSA_Graphs_{objName}")
        tfg.save_obj(DSA_Mean, f"./outputObj/DSA_mean_{objName}")

    # print([legRetList, DSA_Tables, DSA_Graphs, DSA_Mean])
    # del result
    r("rm(list = ls())")
    r("gc()")
    return legRetList, DSA_Tables, DSA_Graphs, DSA_Mean


def deCirces(obj, toName, fromName, legRet, DSA_Table, DSA_Graph, DE_Recp, path, threshold=0.1):
    r["source"]("Codes/Ligand_Receptor_pipeline.R")
    findMarkers_python = r("findMarkers_python")
    DElegenedNcolor = r("DElegenedNcolor")
    DSAlegenedNcolor = r("DSAlegenedNcolor")
    createCircosPlots = r("createCircosPlots")

    with localconverter(default_converter + rpyp.converter):
        markerallL = findMarkers_python(obj, fromName, toName, genes_to_use=list(legRet["Ligand"]), threshold=threshold)
        markerallR = findMarkers_python(obj, toName, fromName, genes_to_use=list(legRet["Receptor"]),
                                        threshold=threshold)

        markerallL["isL"] = list(map(lambda x: x in list(legRet["Ligand"]), markerallL.index))
        markerallR["isR"] = list(map(lambda x: x in list(legRet["Receptor"]), markerallR.index))

        markerallL = markerallL.loc[markerallL["isL"], :]
        markerallR = markerallR.loc[markerallR["isR"], :]

        if len(markerallL.index) > 0 and len(markerallR.index) > 0:
            DE = DElegenedNcolor(obj, fromName, toName, legRet, markerallL, markerallR)
            DSA = DSAlegenedNcolor(DSA_Table)
            toExpression, fromExpression = ExpTables(obj, toName, fromName)

            createCircosPlots(toExpression, legRet, DE, ro.ListVector(DSA_Graph.Tfs_Per_Recp_in_max_flow()), DSA,
                              fromName, toName, "de", obj,
                              DE_Recp, path=path)


def RUN_DE_LR(conf, lst_Tr, lst_Co, pathTr, pathCo, max_workers=20):
    r["source"]("Codes/Ligand_Receptor_pipeline.R")
    toVec = list(conf.loc[conf["Var"] == "toVecForDE", "Value"])[0].split()
    fromVec = list(conf.loc[conf["Var"] == "fromVecForDE", "Value"])[0].split()

    objTr = r(f"objTr = readRDS('{pathTr}')")
    objCo = r(f"objCo = readRDS('{pathCo}')")

    toVec = utils.intersection(toVec, lst_Tr[0].keys())
    toVec = utils.intersection(toVec, lst_Co[0].keys())

    name_obj_clf = list(conf.loc[conf["Var"] == "nameObjCLF", "Value"])[0]
    plotpathT = list(conf.loc[conf["Var"] == "DEplotpathT", "Value"])[0]
    plotpathC = list(conf.loc[conf["Var"] == "DEplotpathC", "Value"])[0]

    args = []
    for toName in toVec:
        fromVecTemp = utils.intersection(fromVec, lst_Tr[0][toName].keys())
        fromVecTemp = utils.intersection(fromVecTemp, lst_Co[0][toName].keys())
        for fromName in fromVecTemp:
            if toName != fromName:
                toExpressionTr, fromExpressionTr = ExpTables(objTr, toName, fromName)
                toExpressionCo, fromExpressionCo = ExpTables(objCo, toName, fromName)
                args.append((toExpressionTr, toExpressionCo, fromExpressionTr, fromExpressionCo,
                             f"plotOut/ROC/RocResult{fromName}_{toName}", lst_Tr[0][toName][fromName],
                             lst_Co[0][toName][fromName]))

    # with BoundedProcessPoolExecutor(max_workers=max_workers) as executor:
    with concurrent.futures.ProcessPoolExecutor(max_workers=max_workers) as executor:

        results = list(executor.map(_helper_Classfier, args))

    counter = 0
    up = {}
    down = {}

    for toName in toVec:
        up[toName] = {}
        down[toName] = {}
        for fromName in fromVec:
            if toName != fromName and len(results) > counter:
                if len(results[counter]) > 1:
                    up[toName][fromName] = pd.DataFrame(results[counter][1])
                    up[toName][fromName]["importances value"] = up[toName][fromName]["importances value"] / \
                                                                  up[toName][fromName]["importances value"].max()
                    up[toName][fromName]["importances value"] = up[toName][fromName]["importances value"].apply(
                        lambda x: x if x > 0.3 else 0.3)

                if len(results[counter]) > 2:
                    down[toName][fromName] = pd.DataFrame(results[counter][2])
                    down[toName][fromName]["importances value"] = down[toName][fromName]["importances value"] / \
                                                                  down[toName][fromName]["importances value"].max()
                    down[toName][fromName]["importances value"] = down[toName][fromName]["importances value"].apply(lambda x: x if x > 0.3 else 0.3)
            counter += 1

    tfg.save_obj(up, f"outputObj/up_{name_obj_clf}")
    tfg.save_obj(down, f"outputObj/down_{name_obj_clf}")

    # with BoundedProcessPoolExecutor(max_workers=max_workers) as executor:
    # with concurrent.futures.ProcessPoolExecutor(max_workers=max_workers) as executor:

    for toName in toVec:
        fromVecTemp = utils.intersection(fromVec, down[toName].keys())
        fromVecTemp = utils.intersection(fromVecTemp, up[toName].keys())
        for fromName in fromVecTemp:
            if fromName != toName:
                if len(down[toName][fromName].index) > 0:
                    deCirces(objCo, toName, fromName, lst_Co[0][toName][fromName],
                             lst_Co[1][toName][fromName],
                             lst_Co[2].return_GD_from_key(toName).return_GD_from_key(fromName),
                             down[toName][fromName], plotpathC)

                if len(up[toName][fromName].index) > 0:
                    # executor.submit
                    deCirces(objTr, toName, fromName, lst_Tr[0][toName][fromName],
                             lst_Tr[1][toName][fromName],
                             lst_Tr[2].return_GD_from_key(toName).return_GD_from_key(fromName),
                             up[toName][fromName], plotpathT)

    return up, down


def build_or_use_classfier(toExp1, toExp2, fromExp1, fromExp2, RocPlotName,
                           lr1, lr2=None, modle=None):
    r["source"]("Codes/Ligand_Receptor_pipeline.R")
    ProtNamesInfo = pd.read_csv("files/humanProtinInfo.csv")
    ProtActions = pd.read_csv("files/humanProtinAction.csv")
    lst = pyd.main_py_to_R(ProtNamesInfo, ProtActions)
    if lr2 is not None:
        legRet = pd.concat([lr1, lr2])
    else:
        legRet = lr1
    with localconverter(default_converter + rpyp.converter):
        DSA_Table_1 = tfg.DSA_anaylsis(toExp1, legRet["Receptor"], ProtActions, ProtNamesInfo, tfs=lst[0])
        DSA_Table_2 = tfg.DSA_anaylsis(toExp2, legRet["Receptor"], ProtActions, ProtNamesInfo, tfs=lst[0])

    args1 = [toExp1, fromExp1, DSA_Table_1[0], legRet]
    args2 = [toExp2, fromExp2, DSA_Table_2[0], legRet]

    if modle is None:
        print("Stop0")
        return utils.DSA_Classfier(args1, args2, return_modle=False, plot_name=RocPlotName)
    return utils.DSA_Classfier(args1, args2, modle=modle, return_modle=False, plot_name=RocPlotName)


def main():
    try:
        conf = pd.read_csv("config.csv")
    except:
        raise Exception("no input CSV found")

    pathTret = list(conf.loc[conf["Var"] == "objTr", "Value"])[0]
    pathControl = list(conf.loc[conf["Var"] == "objCo", "Value"])[0]

    if type(pathControl) == float and np.isnan(pathControl):
        flag = False
    else:
        flag = True

    plotpath = list(conf.loc[conf["Var"] == "plotpathTr", "Value"])[0]
    toVec = list(conf.loc[conf["Var"] == "toVec", "Value"])[0].split(" ")
    fromVec = list(conf.loc[conf["Var"] == "fromVec", "Value"])[0].split(" ")
    args = [(f"InputObj/{pathTret}", toName, fromName, plotpath) for toName in toVec for fromName in fromVec if
            fromName != toName]

    name_obj_Tr = list(conf.loc[conf["Var"] == "nameObjTr", "Value"])[0]
    lst_Tr = run_pipline(args, name_obj_Tr)
    if flag:
        plotpath = list(conf.loc[conf["Var"] == "plotpathCo", "Value"])[0]
        args = [(f"InputObj/{pathControl}", toName, fromName, plotpath) for toName in toVec for fromName in fromVec if
                fromName != toName]

        name_obj_Co = list(conf.loc[conf["Var"] == "nameObjCo", "Value"])[0]
        lst_Co = run_pipline(args, name_obj_Co)
        up, down = RUN_DE_LR(conf, lst_Tr, lst_Co, f"InputObj/{pathTret}", f"InputObj/{pathControl}")

        projName = list(conf.loc[conf["Var"] == "Proj", "Value"])[0]
        graph = neo.start_connection()
        neo.delete_nodes_by_prop(graph, prop="project", key=projName, delete_all_grpah=False)

        for toName in toVec:
            toNodeTr = neo.add_node_to_neo(graph, toName, name_obj_Tr, projName)
            toNodeCo = neo.add_node_to_neo(graph, toName, name_obj_Co, projName)

            neo.crate_all_from_nodes(graph, toNodeTr, lst_Tr[3][toName],
                                     lst_Tr[1][toName], lst_Tr[0][toName], lst_Tr[2].return_GD_from_key(toName),
                                     name_obj_Tr, projName)

            neo.crate_all_from_nodes(graph, toNodeCo, lst_Co[3][toName],
                                     lst_Co[1][toName], lst_Co[0][toName], lst_Co[2].return_GD_from_key(toName),
                                     name_obj_Co, projName)

            neo.make_cluster_trat_relationship(graph, toNodeTr, toNodeCo, down[toName])

            neo.make_cluster_trat_relationship(graph, toNodeCo, toNodeTr, up[toName])

        return lst_Tr, lst_Co

    else:

        projName = list(conf.loc[conf["Var"] == "Proj", "Value"])[0]
        graph = neo.start_connection()
        neo.delete_nodes_by_prop(graph, prop="project", key=projName, delete_all_grpah=False)

        for toName in toVec:
            toNodeTr = neo.add_node_to_neo(graph, toName, name_obj_Tr, projName)
            neo.crate_all_from_nodes(graph, toNodeTr, lst_Tr[3][toName],
                                     lst_Tr[1][toName], lst_Tr[0][toName], lst_Tr[2].return_GD_from_key(toName),
                                     name_obj_Tr, projName)

    return lst_Tr


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


def DE_DSA():
    conf = pd.read_csv("config.csv")

    name_obj_Tr = list(conf.loc[conf["Var"] == "nameObjTr", "Value"])[0]
    name_obj_Co = list(conf.loc[conf["Var"] == "nameObjCo", "Value"])[0]

    DSA_MeanTr = tfg.load_obj(f"outputObj/DSA_mean_{name_obj_Tr}")
    DSA_MeanCo = tfg.load_obj(f"outputObj/DSA_mean_{name_obj_Co}")

    DSA_UP = pd.DataFrame()
    DSA_DOWN = pd.DataFrame()

    for key, value in DSA_MeanTr.items():
        value["To Cluster"] = key
        value = value[["Cluster", "To Cluster", "DSA_Mean"]]
        value.columns = (["From Cluster", "To Cluster", "DSA_Mean"])
        DSA_UP = pd.concat([DSA_UP, value])

    for key, value in DSA_MeanCo.items():
        value["To Cluster"] = key
        value = value[["Cluster", "To Cluster", "DSA_Mean"]]
        value.columns = (["From Cluster", "To Cluster", "DSA_Mean"])
        DSA_DOWN = pd.concat([DSA_DOWN, value])

    DSA_UP["From Cluster"] = DSA_UP["From Cluster"].apply(format_dsa_file)
    DSA_UP["To Cluster"] = DSA_UP["To Cluster"].apply(format_dsa_file)

    DSA_DOWN["From Cluster"] = DSA_DOWN["From Cluster"].apply(format_dsa_file)
    DSA_DOWN["To Cluster"] = DSA_DOWN["To Cluster"].apply(format_dsa_file)

    DSA_DOWN.to_csv(r"./files/DSA_DOWN.csv")
    DSA_UP.to_csv(r"./files/DSA_UP.csv")

    return DSA_UP, DSA_DOWN


def test_function():
    conf = pd.read_csv(("config.csv"))
    pathTret = list(conf.loc[conf["Var"] == "objTr", "Value"])[0]
    pathControl = list(conf.loc[conf["Var"] == "objCo", "Value"])[0]
    toVec = list(conf.loc[conf["Var"] == "toVec", "Value"])[0].split(" ")
    fromVec = list(conf.loc[conf["Var"] == "fromVec", "Value"])[0].split(" ")
    name_obj_Tr = list(conf.loc[conf["Var"] == "nameObjTr", "Value"])[0]
    name_obj_Co = list(conf.loc[conf["Var"] == "nameObjCo", "Value"])[0]
    name_obj_clf = list(conf.loc[conf["Var"] == "nameObjCLF", "Value"])[0]
    plotpathT = list(conf.loc[conf["Var"] == "DEplotpathT", "Value"])[0]
    plotpathC = list(conf.loc[conf["Var"] == "DEplotpathC", "Value"])[0]

    legRetListTr = tfg.load_obj(f"outputObj/legRetLists_{name_obj_Tr}")
    DSA_TablesTr = tfg.load_obj(f"outputObj/DSA_Tables_{name_obj_Tr}")
    DSA_GraphsTr = tfg.load_obj(f"outputObj/DSA_Graphs_{name_obj_Tr}")
    DSA_MeanTr = tfg.load_obj(f"outputObj/DSA_mean_{name_obj_Tr}")
    lst_Tr = (legRetListTr, DSA_TablesTr, DSA_GraphsTr, DSA_MeanTr)

    legRetListCo = tfg.load_obj(f"outputObj/legRetLists_{name_obj_Co}")
    DSA_TablesCo = tfg.load_obj(f"outputObj/DSA_Tables_{name_obj_Co}")
    DSA_GraphsCo = tfg.load_obj(f"outputObj/DSA_Graphs_{name_obj_Co}")
    DSA_MeanCo = tfg.load_obj(f"outputObj/DSA_mean_{name_obj_Co}")
    lst_Co = (legRetListCo, DSA_TablesCo, DSA_GraphsCo, DSA_MeanCo)

    up = tfg.load_obj(f"outputObj/up_{name_obj_clf}")
    down = tfg.load_obj(f"outputObj/up_{name_obj_clf}")

    objTr = r(f"objTr = readRDS('InputObj/{pathTret}')")
    objCo = r(f"objCo = readRDS('InputObj/{pathControl}')")

    # toVec = intersection(toVec, lst_Tr[0].keys())
    # toVec = intersection(toVec, lst_Co[0].keys())
    # projName = list(conf.loc[conf["Var"] == "Proj", "Value"])[0].split()

    for toName in toVec:
        fromVecTemp = utils.intersection(fromVec, lst_Tr[0][toName].keys())
        fromVecTemp = utils.intersection(fromVecTemp, lst_Tr[0][toName].keys())
        for fromName in fromVecTemp:
            if fromName != toName:
                if len(down[toName][fromName].index) > 0:
                    deCirces(objCo, toName, fromName, lst_Co[0][toName][fromName],
                             lst_Co[1][toName][fromName],
                             lst_Co[2].return_GD_from_key(toName).return_GD_from_key(fromName),
                             list(down[toName][fromName]["feture"]), plotpathC)

                if len(up[toName][fromName].index) > 0:
                    # executor.submit
                    deCirces(objTr, toName, fromName, lst_Tr[0][toName][fromName],
                             lst_Tr[1][toName][fromName],
                             lst_Tr[2].return_GD_from_key(toName).return_GD_from_key(fromName),
                             list(up[toName][fromName]["feture"]), plotpathT)


if __name__ == "__main__":
    print(os.getcwd())
    main()
# test_function()

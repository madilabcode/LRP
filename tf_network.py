import numpy as np
import pandas as pd
import os
os.environ["R_HOME"] = r"C:\Program Files\R\R-4.1.1"
os.environ['path'] += r";C:\Program Files\R\R-4.1.1\bin"
from rpy2.robjects import r
from rpy2.robjects.conversion import localconverter
from rpy2.robjects import default_converter
from Codes import Tf_Graph as tfg
from Codes import utils
import rpy2.robjects.pandas2ri as rpyp
import concurrent.futures

def tfs_importance(tr=True):
    try:
        conf = pd.read_csv("config.csv")

    except:
        raise Exception("no input CSV found")

    name_obj_Tr = list(conf.loc[conf["Var"] == "nameObjTr", "Value"])[0]
    name_obj_Co = list(conf.loc[conf["Var"] == "nameObjCo", "Value"])[0]

    if tr:
        graphs = tfg.load_obj(f"outputObj\DSA_Graphs_{name_obj_Tr}")
    else:
        graphs = tfg.load_obj(f"outputObj\DSA_Graphs_{name_obj_Co}")
    toCluster = input("To Cluster")
    fromCluster = input("From Cluster")
    graphObj = graphs.return_GD_from_key(toCluster).return_GD_from_key(fromCluster).update_obj()
    df = graphObj.calculate_significant_tf_Multy_sinc()
    df.to_csv(f"./files/{fromCluster}_{toCluster}_{tr}_importance.csv")


def tf_flow_graph(tr=True):
    r["source"]('Codes/Ligand_Receptor_pipeline.R')
    try:
        conf = pd.read_csv("config.csv")

    except:
        raise Exception("no input CSV found")

    name_obj_Tr = list(conf.loc[conf["Var"] == "nameObjTr", "Value"])[0]
    name_obj_Co = list(conf.loc[conf["Var"] == "nameObjCo", "Value"])[0]

    if tr:
        graphs = tfg.load_obj(f"outputObj\DSA_Graphs_{name_obj_Tr}")
    else:
        graphs = tfg.load_obj(f"outputObj\DSA_Graphs_{name_obj_Co}")

    toCluster = input("To Cluster")
    fromCluster = input("From Cluster")
    tfs = input("enter all tfs")
    tfsList = tfs.split()
    graphObj = graphs.return_GD_from_key(toCluster).return_GD_from_key(fromCluster).update_obj()

    for tf in tfsList:
        try:
            df = graphObj.flow_dict_to_single_tf(tf)
            df["flow"] = df["flow"].apply(lambda x: round(x, 3))
            df["flow"] = df["flow"].astype(float)
            df = df.loc[df.flow > 0.1]
            df["flow"] = df["flow"].astype(str)
            df.to_csv(f"./files/{toCluster}_{fromCluster}_{tf}_flow.csv")
              #  draw_flow_graph_from_df(df, root=tf)
                    #input("Press Enter to continue...")


        except:
            print("tf not in grpah")


if __name__ == "__main__":
    func = input("You want to make tf flow graph [f] or calc tf importance [i]: ")
    tr = input("is tr obj? [y/n]")
    if func == "i":
        if tr == "y":
            tfs_importance()
        else:
            tfs_importance(False)
    else:
        if tr == "y":
          tf_flow_graph()
        else:
            tf_flow_graph(False)

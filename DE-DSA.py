import numpy as np
import pandas as pd
from bounded_pool_executor import BoundedProcessPoolExecutor
#import anndata2ri as ri
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


def format_dsa_file(x):
    annot = pd.read_csv("Annot.csv")
    try:
      type_cluster = list(annot.loc[annot['Cluster'] == int(x),'Type'])[0]
      subtype_cluster = list(annot.loc[annot['Cluster'] == int(x),'subType'])[0]
    except :
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
        value.columns =(["From Cluster","To Cluster", "DSA_Mean"])
        DSA_UP = pd.concat([DSA_UP, value])

    for key, value in DSA_MeanCo.items():
        value["To Cluster"] = key
        value = value[["Cluster", "To Cluster", "DSA_Mean"]]
        value.columns = (["From Cluster","To Cluster", "DSA_Mean"])
        DSA_DOWN = pd.concat([DSA_DOWN, value])

    DSA_UP["From Cluster"] = DSA_UP["From Cluster"].apply(format_dsa_file)
    DSA_UP["To Cluster"] = DSA_UP["To Cluster"].apply(format_dsa_file)

    DSA_DOWN["From Cluster"] = DSA_DOWN["From Cluster"].apply(format_dsa_file)
    DSA_DOWN["To Cluster"] = DSA_DOWN["To Cluster"].apply(format_dsa_file)

    DSA_DOWN.to_csv(r"./files/DSA_DOWN.csv")
    DSA_UP.to_csv(r"./files/DSA_UP.csv")

    return DSA_UP, DSA_DOWN
    
    
if __name__ == "__main__":
  DE_DSA()

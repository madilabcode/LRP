import pandas as pd
import numpy as np
import random
import math
import scipy.stats as st
from Codes import randomForstSIngleClassfier as rf
import zlib, json, base64
from rpy2.robjects import r
from rpy2 import robjects as ro
from rpy2.robjects.conversion import localconverter
from rpy2.robjects import default_converter
import rpy2.robjects.pandas2ri as rpyp


def intersection(lst1, lst2):
    return [value for value in lst1 if value in lst2]


def align_clusters(orig, active1, active2, name1, name2, format1=None, format2=None):
    format1 = dict(format1)
    format2 = dict(format2)
    orig["cell"] = orig.index
    active1["cell"] = active1.index
    active2["cell"] = active2.index

    if format1 is not None:
        active1["cell"] = active1["cell"].apply(
            lambda cell: format1[cell.split("_")[0]] + cell.split("_")[len(cell.split("_")) - 1])

    if format2 is not None:
        active2["cell"] = active2["cell"].apply(
            lambda cell: format2[cell.split("_")[0]] + cell.split("_")[len(cell.split("_")) - 1])

    orig.index = range(len(orig.index))
    active1.index = range(len(active1.index))
    active2.index = range(len(active2.index))

    orig = orig.rename(columns={orig.columns[0]: 'orig_ident'})
    active1 = active1.rename(columns={active1.columns[0]: 'active1_ident'})
    active2 = active2.rename(columns={active2.columns[0]: 'active2_ident'})

    active1 = pd.merge(orig, active1, how="inner", on="cell")
    active2 = pd.merge(orig, active2, how="inner", on="cell")

    active1_result = pd.pivot_table(active1, index=["active1_ident", "orig_ident"], values=["cell"], aggfunc=len)
    active2_result = pd.pivot_table(active2, index=["active2_ident", "orig_ident"], values=["cell"], aggfunc=len)

    active1_result["active1_ident"] = list(map(lambda ind: ind[0], active1_result.index))
    active1_result["orig_ident"] = list(map(lambda ind: ind[1], active1_result.index))
    active1_result = active1_result[["active1_ident", "orig_ident", "cell"]]
    active1_result.index = range(len(active1_result.index))

    active2_result["active2_ident"] = list(map(lambda ind: ind[0], active2_result.index))
    active2_result["orig_ident"] = list(map(lambda ind: ind[1], active2_result.index))
    active2_result = active2_result[["active2_ident", "orig_ident", "cell"]]
    active2_result.index = range(len(active2_result.index))

    active1_sum_of_cluster = pd.pivot_table(active1, index=["active1_ident"], values=["cell"], aggfunc=len)
    active1_sum_of_cluster["active1_ident"] = active1_sum_of_cluster.index
    active1_sum_of_cluster.index = range(len(active1_sum_of_cluster.index))
    active2_sum_of_cluster = pd.pivot_table(active2, index=["active2_ident"], values=["cell"], aggfunc=len)
    active2_sum_of_cluster["active2_ident"] = active2_sum_of_cluster.index
    active2_sum_of_cluster.index = range(len(active2_sum_of_cluster.index))

    active1_result = pd.merge(active1_result, active1_sum_of_cluster, on="active1_ident")
    active1_result["P_of_cluster"] = active1_result["cell_x"] / active1_result["cell_y"]
    active1_result = active1_result.loc[active1_result["P_of_cluster"] > 0.5, :].reindex(
        columns=["active1_ident", "orig_ident", "P_of_cluster"])
    active2_result = pd.merge(active2_result, active2_sum_of_cluster, on="active2_ident")
    active2_result["P_of_cluster"] = active2_result["cell_x"] / active2_result["cell_y"]
    active2_result = active2_result.loc[active2_result["P_of_cluster"] > 0.5, :].reindex(
        columns=["active2_ident", "orig_ident", "P_of_cluster"])

    result = pd.merge(active1_result, active2_result, on="orig_ident").reindex(
        columns=["active1_ident", "active2_ident"])
    result = result.rename(columns={"active1_ident": name1, "active2_ident": name2})
    return result


def reduce_table_by_key(table, key, Col_Name):
    table["is_Key"] = table[Col_Name].apply(lambda name: 1 if str(name).find(key) != -1 else 0)
    table = table.loc[table["is_Key"] == 1, :]
    table = table.drop(axis=1, labels="is_Key")
    return table


def sum_avg_ligand(ligs, exp):
    sub_exp = exp.loc[exp.index.isin(ligs), :]
    sub_exp_Avg = sub_exp.apply(np.mean, axis=1)
    return np.log1p(np.sum(sub_exp_Avg))


def dsa_table_update_to_lig(toExp, fromExp, DSA_Table, legRet):
    toExp = toExp.loc[toExp.index.isin(np.unique(DSA_Table["Recp"])), :]

    for rec in DSA_Table.index:
        sum_lig = sum_avg_ligand(list(legRet.loc[legRet["Receptor"] == rec, "Ligand"]), fromExp)
        DSA_Table.loc[rec, "DSA"] *= sum_lig

    for gene in toExp.index:
        toExp.loc[gene, :] = np.log1p(toExp.loc[gene, :]) * DSA_Table.loc[gene, "DSA"]

    return toExp


def dsa_with_lig(toExp, fromExp, DSA_Table, legRet):
    toExp = toExp.loc[toExp.index.isin(np.unique(DSA_Table["Recp"])), :]

    for rec in DSA_Table["Recp"]:
        sum_lig = sum_avg_ligand(list(legRet.loc[legRet["Receptor"] == rec, "Ligand"]), fromExp)
        DSA_Table.loc[DSA_Table["Recp"] == rec, "DSA"] *= sum_lig

    for gene in toExp.index:
        toExp.loc[gene, :] = np.log1p(toExp.loc[gene, :]) * float(DSA_Table.loc[DSA_Table["Recp"] == gene, "DSA"])

    return toExp.apply(np.sum, axis=0)


def DSA_DE(dsa1, dsa2):
    # dsa1 = dsa_with_lig(args1)

    dsa1 = list(filter(lambda x: not math.isnan(x), list(dsa1)))
    dsa2 = list(filter(lambda x: not math.isnan(x), list(dsa2)))
    stat, pval = st.ttest_ind(dsa1, dsa2)
    return stat, pval


def DSA_Classfier(args1, args2, modle=None, return_modle=False, plot_name="result"):
    args1 = tuple(args1)
    args2 = tuple(args2)

    dsa_t1 = dsa_table_update_to_lig(*args1)
    dsa_t2 = dsa_table_update_to_lig(*args2)

    dsa_t1 = dsa_t1.transpose()
    dsa_t2 = dsa_t2.transpose()

    dsa_t1["ident"] = 1
    dsa_t2["ident"] = 0

    if len(dsa_t1.index) > len(dsa_t2.index):
        rows = random.sample(list(dsa_t1.index), len(dsa_t2.index))
        dsa_t1 = dsa_t1.loc[dsa_t1.index.isin(rows), :]
    else:
        rows = random.sample(list(dsa_t2.index), len(dsa_t1.index))
        dsa_t2 = dsa_t2.loc[dsa_t2.index.isin(rows), :]

    dsa_t1.index = range(len(dsa_t1.index))
    dsa_t2.index = range(len(dsa_t2.index))

    dsa_table = pd.concat([dsa_t1, dsa_t2], axis=0)

    for col in dsa_table.columns:
        dsa_table[col] = dsa_table[col].apply(lambda x: 0 if math.isnan(x) else x)

    dsa_table.index = range(len(dsa_table.index))

    if modle is not None:
        ls = rf.use_ex_modle(modle, dsa_table, plot_name)
    else:
        ls = rf.random_forst_expresstion(dsa_table, plot_name=plot_name, returnModle=return_modle)

    if return_modle:
        return dsa_table, ls[0], ls[1], ls[2]

    return dsa_table, ls[0], ls[1]


def shift_dist_of_object(obj_tr, obj_co, assay="scale.data"):
    with localconverter(default_converter + rpyp.converter):
        r("library(Seurat)")
        r("library(dplyr)")
        GetAssayData = r(f"function(obj) obj[['RNA']]@{assay} %>% as.data.frame()")
        shift_exp = r("""function(obj,exp){ obj[['RNA']]@scale.data =  as.matrix(exp)
         return(obj)}""")
        exp_tr, exp_co = GetAssayData(obj_tr), GetAssayData(obj_co)
        min_value_tr = exp_tr.apply(min, axis=1).min()
        min_value_co = exp_co.apply(min, axis=1).min()
        min_value = min(min_value_tr, min_value_co)
        exp_tr += np.abs(min_value)
        exp_co += np.abs(min_value)
        return shift_exp(obj_tr, exp_tr), shift_exp(obj_co, exp_co)


def save_to_csv(table, name):
    table.to_csv(name)


def json_zip(js, ZIPJSON_KEY='base64(zip(o))'):
    js = {
        ZIPJSON_KEY: base64.b64encode(
            zlib.compress(
                json.dumps(js).encode('utf-8')
            )
        ).decode('ascii')
    }
    return json.dumps(js)


def json_unzip(zp, ZIPJSON_KEY='base64(zip(o))'):
    zp = json.loads(zp)
    assert (ZIPJSON_KEY in zp.keys())
    try:
        js = zlib.decompress(base64.b64decode(zp[ZIPJSON_KEY]))
    except:
        raise RuntimeError("Could not decode/unzip the contents")
    try:
        js = json.loads(js)
    except:
        raise RuntimeError("Could interpret the unzipped contents")
    return js

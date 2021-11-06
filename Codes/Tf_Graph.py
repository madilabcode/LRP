import pandas as pd
import numpy as np
import math
from rpy2.robjects import r
# import matplotlib.pyplot as plt
from sklearn.metrics import roc_curve, roc_auc_score
from sklearn.ensemble import RandomForestClassifier
import scipy.stats as sc
import networkx as nx
import pickle
from Codes import CERNO as ce
import concurrent.futures

CUTOFF = 0.1
class Obj_dict:
    def __init__(self):
        self.dic = {}

    def add_obj(self, key, obj):
        self.dic[key] = obj

    def return_GD_from_key(self, key):
        return self.dic[key]

    def save_obj(self, name):
        save_obj(self, name)


class graphs_dict:
    def __init__(self, flow_dics, capcity_network):
        self.flow_dics = flow_dics
        self.capcity_network = capcity_network
        self.max_flow_dict = self.calculate_max_flow_for_all_recp()
        self.max_multy_flow = None

    def sub_capcity_gaph_to_dict(self):
        sub = nx.to_pandas_edgelist(self.capcity_network)
        sub = sub.loc[sub["capacity"] < np.inf, :]
        gp = nx.from_pandas_edgelist(sub, "source", "target", edge_attr="capacity", create_using=nx.DiGraph())
        return nx.to_dict_of_dicts(gp)

    @classmethod
    def make_capcity_graph(cls, network_dict):
        ProtInfo = pd.read_csv(r"files/humanProtinInfo.csv")
        pa = pd.read_csv(r"files/humanProtinAction.csv")
        pa["Target"] = pa["Output-node Gene Symbol"]
        pa["Source"] = pa["Input-node Gene Symbol"]
        pa["capacity"] = math.inf
        keeps = ["Source", "Target", "capacity"]
        pa = pa[keeps]
        genes = network_dict.keys()
        genes = list(map(lambda gene: gene.split("_")[0], genes))

        pa["Source"] = pa["Source"].apply(format_Gene_dict).astype(str)
        pa["Target"] = pa["Target"].apply(format_Gene_dict).astype(str)

        pa = pa.loc[pa["Source"].isin(genes), :]
        pa = pa.loc[pa["Target"].isin(genes), :]
        pa["Source"] += "_Source"
        pa["Target"] += "_Target"

        gp1 = nx.from_pandas_edgelist(pa, "Source", "Target", edge_attr="capacity", create_using=nx.DiGraph())
        gp2 = nx.from_dict_of_dicts(network_dict, create_using=nx.DiGraph())

        return nx.compose(gp1, gp2)

    @classmethod
    def make_obj_out_of_capcity_network(cls, network_dict, recptors):
        gp = graphs_dict.make_capcity_graph(network_dict)
        flow_dicts = {}

        recptors = np.unique(recptors)
        for rec in recptors:
            try:
                flow_value, flow_dict = nx.maximum_flow(gp, "s", rec + "_Source", flow_func=nx.algorithms.flow.dinitz)
                flow_dicts[rec] = flow_dict
            except:
                print("node " + rec + " is not in the graph")

        return graphs_dict(flow_dicts, gp)

    def update_obj(self):
        return graphs_dict(self.flow_dics, self.capcity_network)

    def calculate_max_flow(self, recp):
        dfFlow = self.make_df_flow(recp)
        max_flow = 0

        for row in dfFlow.iterrows():
            row = row[1]
            if row["target"] == "sink" :
                max_flow += row["flow"]

        return max_flow

    def calculate_max_flow_for_all_recp(self):
        max_flow_dict = {}
        for recp in self.flow_dics.keys():
            max_flow_dict[recp] = self.calculate_max_flow(recp)

        return max_flow_dict

    def make_df_flow(self, recp=None, flow=None):
        if flow is None:
            flow = self.flow_dics[recp]
        source = []
        target = []
        flows = []

        for Inode in flow.keys():
            edges = flow[Inode]
            for Onode in edges.keys():
                if edges[Onode] != 0 and Inode != Onode:
                    source.append(Inode)
                    target.append(Onode)
                    flows.append(edges[Onode])

        df = pd.DataFrame([source, target, flows])
        df = df.transpose()
        df.columns = ["source", "target", "flow"]
        return df

    def Tfs_Per_Recp_in_max_flow(self):
        tf_rec = {}
        for rec in self.flow_dics.keys():
            df = self.make_df_flow(rec)
            df = df.loc[df["target"] == "sink", :]
            tf_rec[rec] = len(df["source"])

        return tf_rec

    def calculate_significant_tf(self, recp):
        df = self.make_df_flow(recp)
        all_tf = np.unique(df.loc[df["target"] == "sink", "source"])
        tf_effects = {"tfs": [], "flow_delta": []}

        for tf in all_tf:
            pa = nx.to_pandas_edgelist(self.capcity_network)
            pa.columns = ["source", "target", "capacity"]
            pa["isNotTf"] = pa.apply(
                lambda row: True if str(row["source"]).find(tf) < 0 and str(row["target"]).find(tf) < 0 else False,
                axis=1)
            pa = pa.loc[pa["isNotTf"], :]
            gpf = nx.from_pandas_edgelist(pa, "source", "target", edge_attr="capacity", create_using=nx.DiGraph())
            flow_value = nx.maximum_flow(gpf, recp, "sink", flow_func=nx.algorithms.flow.dinitz)[0]
            tf_effects["tfs"].append(tf)
            tf_effects["flow_delta"].append(self.max_flow_dict[recp] - flow_value)

        tf_effects = pd.DataFrame.from_dict(tf_effects)
        tf_effects = tf_effects.sort_values(by="flow_delta", axis=0, ascending=False)
        return tf_effects

    @staticmethod
    def _calculate_flow_delta_for_df(args):
        gene, pat, max_flow = args
        if gene != "s" and gene != "t":
            pat["isNotTf"] = pat.apply(
                lambda row: True if str(row["source"]).find(gene) < 0 and str(row["target"]).find(
                    gene) < 0 else False, axis=1)
            pat = pat.loc[pat["isNotTf"], :]
            gpf = nx.from_pandas_edgelist(pat, "source", "target", edge_attr="capacity", create_using=nx.DiGraph())
            flow_value = nx.maximum_flow(gpf, "s", "t", flow_func=nx.algorithms.flow.dinitz)[0]
            return gene, max_flow - flow_value

    def calculate_significant_tf_Multy_sinc(self, only_tf=True):
        pa = nx.to_pandas_edgelist(self.capcity_network)
        pa.columns = ["source", "target", "capacity"]
        tf_effects = {"tfs": [], "flow_delta": []}
        all_genes = []
        n = len(pa.index)

        for recp in self.flow_dics.keys():
            df = self.make_df_flow(recp)
            if only_tf:
                all_genes += list(df.loc[df["target"] == "s", "source"])
            else:
                all_genes += list(df["source"]) + list(df["target"])
            row = [recp + "_Source", "t", math.inf]
            pa.loc[n] = row
            n += 1

        all_genes = np.unique(all_genes)
        gpf = nx.from_pandas_edgelist(pa, "source", "target", edge_attr="capacity", create_using=nx.DiGraph())
        max_flow = nx.maximum_flow(gpf, "s", "t", flow_func=nx.algorithms.flow.dinitz)[0]
        self.max_multy_flow = max_flow

        with concurrent.futures.ProcessPoolExecutor(max_workers=4) as executor:
            results = list(executor.map(graphs_dict._calculate_flow_delta_for_df,
                                        [(gene, pa.copy(), max_flow) for gene in all_genes]))

        for element in results:
            tf_effects["tfs"].append(element[0])
            tf_effects["flow_delta"].append(element[1])

        tf_effects = pd.DataFrame.from_dict(tf_effects)
        tf_effects = tf_effects.sort_values(by="flow_delta", axis=0, ascending=False)
        return tf_effects

    def run_multiy_source_flow(self,graph):
        n = len(graph.index)
        for recp in self.flow_dics.keys():
            df = self.make_df_flow(recp)
            row = ["source_node", recp,np.inf]
            graph.loc[n] = row
            n += 1

        gpf = nx.from_pandas_edgelist(graph, "source", "target", edge_attr="capacity", create_using=nx.DiGraph())
        max_flow, flow = nx.maximum_flow(gpf, "source_node", "sink", flow_func=nx.algorithms.flow.dinitz)

        return  max_flow, flow

    def add_t_to_network(self):
        pa = nx.to_pandas_edgelist(self.capcity_network)
        pa.columns = ["source", "target", "capacity"]
        n = len(pa.index)
        for recp in self.flow_dics.keys():
            row = [recp + "_Source", "t", math.inf]
            pa.loc[n] = row
            n += 1
        return pa


    def calculate_max_flow_for_one_tf(self, tf):
        pa = self.add_t_to_network()
        cols = ["source", "target", "capacity"]
        pa["istf"] = pa.apply(lambda row: True if row["source"] != "s" or row["target"].find(tf) >= 0 else False,
                              axis=1)
        pa = pa.loc[pa["istf"], cols]
        gpf = nx.from_pandas_edgelist(pa, "source", "target", edge_attr="capacity", create_using=nx.DiGraph())
        max_flow, flow = self.run_multiy_source_flow(pa)
        result = self.make_df_flow(flow=flow)
        result["notTorS"] = result.apply(lambda row: True if row["source"] != "source_node" and row["target"] != "sink" else False,
                                         axis=1)
        return result.loc[result["notTorS"], ["source", "target", "flow"]]

    def calculate_p_value_for_recp(self, rec, num_of_perm=1000, fix_rec=True):
        flow_list = []
        stat = nx.maximum_flow(self.capcity_network, "s", rec + "_Source", flow_func=nx.algorithms.flow.dinitz)[0]
        graph = nx.to_pandas_edgelist(self.capcity_network)


def shortest_path(grpah, node1, node2, upload_name="temp_grpah"):
    try:
        return nx.single_source_dijkstra(grpah, node1, node2)
    except:
        grpah = load_obj(upload_name)
        return nx.single_source_dijkstra(grpah, node1, node2)


def format_Gene_dict(gene, from_dict=False):
    if from_dict:
        gene = gene[1:len(gene) - 1]

    fl = gene[0]
    sl = gene[1:]
    sl = str.lower(sl)

    return fl + sl


def find_wights(edge):
    return (1 / edge["inputAvgExp"]) * (1 / edge["OutputAvgExp"]) * (2 - edge["weight"])


def save_obj(obj, name):
    with open(name + '.pkl', 'wb') as f:
        pickle.dump(obj, f, pickle.HIGHEST_PROTOCOL)


def load_obj(name):
    with open(name + '.pkl', 'rb') as f:
        return pickle.load(f)


def hip_Test(targets, N, exp_cell):
    k = len(exp_cell)
    target_exp = exp_cell[exp_cell.index.isin(targets)]
    m = sum(target_exp)
    q = len(target_exp)
    n = N - m
    p_value = float(r(f"python_phyper({q - 1}, {m}, {n}, {k})"))
    return (p_value)


def hiper_test_tfs(exp, tfs):
    tf_score = {}
    clusterExprsstion = exp.apply(np.sum, axis=1)
    clusterExprsstion = clusterExprsstion[clusterExprsstion > 0]
    for tf in tfs.keys():
        try:
            score = hip_Test(tfs[tf], sum(clusterExprsstion), clusterExprsstion)
        except:
            score = hip_Test([tfs[tf]], sum(clusterExprsstion), clusterExprsstion)

        if score < 0.05:
            tf_score[tf] = score

    return tf_score


def Tf_Per_Recp(exp, recptors, tfs):
    recptors = np.unique(recptors)
    recp_tf = {recp: 0 for recp in recptors}
    clusterExprsstion = exp.apply(np.sum, axis=1)
    clusterExprsstion = clusterExprsstion[clusterExprsstion > 0]

    for tf in tfs.keys():
        for targ in tfs[tf]:
            if targ in clusterExprsstion.index and targ in recptors:
                recp_tf[targ] += 1

    return recp_tf


def build_flowing_network_with_normlized_wights(ProtAction, tfs_enrc, exp):
    pa = ProtAction
    pa["Source"] = pa["Input-node Gene Symbol"]
    pa["Target"] = pa["Output-node Gene Symbol"]
    pa["capacity"] = pa["Edge direction score"]
    keeps = ["Source", "Target", "capacity"]
    pa = pa[keeps]

    pa["Source"] = pa["Source"].apply(format_Gene_dict).astype(str)
    pa["Target"] = pa["Target"].apply(format_Gene_dict).astype(str)

    clusterExprsstion = exp.copy()
    clusterExprsstion["not_zero"] = clusterExprsstion.apply(lambda x: x.astype(bool).sum(), axis=1)
    clusterExprsstion = clusterExprsstion.loc[clusterExprsstion.not_zero > clusterExprsstion.shape[1] * CUTOFF, :]
    clusterExprsstion.drop("not_zero",axis=1,inplace = True)
    clusterExprsstion = pd.DataFrame(exp.mean(axis=1))


    pa = pa.loc[pa["Source"].isin(list(clusterExprsstion.index)), :]
    pa = pa.loc[pa["Target"].isin(list(clusterExprsstion.index)), :]
    pa.reset_index(drop=True,inplace=True)

    gene_wights = pd.DataFrame(clusterExprsstion).copy()
    gene_wights["gene"] = clusterExprsstion.index
    gene_wights["exp"] = clusterExprsstion.values
    loc, scale = sc.expon.fit(gene_wights.exp)
    
    gene_wights["wights"] = gene_wights["exp"].apply(lambda x: sc.expon.cdf(x, loc, scale))
    gene_wights = gene_wights[["gene", "wights"]]

    pa["capacity"] = pa.apply(lambda x: float(x["capacity"]) * float(gene_wights.loc[gene_wights.gene == x["Source"], "wights"]) * \
                                        float(gene_wights.loc[gene_wights.gene == x["Target"], "wights"]),axis=1)
    n = pa.shape[0]
    for tf in tfs_enrc:
        if tf in  pa["Target"].to_list():
            row = [tf, "sink", np.inf]
            pa.loc[n] = row
            n += 1

    gp = nx.from_pandas_edgelist(pa, "Source", "Target", edge_attr="capacity", create_using=nx.DiGraph())
    save_obj(gp, "temp_grpah")

    return gp


def DSA_anaylsis(exp, recptors, ProtAction, ProtInfo, tfs):
    # tfs_scores = hiper_test_tfs(exp, tfs)
    tfs_scores = ce.CERNO_alg(exp.apply(np.sum, axis=1), tfs)
    gpf = build_flowing_network_with_normlized_wights(ProtAction, tfs_scores, exp)
    flow_values = []
    flow_dicts = {}
    recptors = np.unique(recptors)

    for rec in recptors:
        try:
            flow_value, flow_dict = nx.maximum_flow(gpf, rec ,"sink", flow_func=nx.algorithms.flow.dinitz)
            flow_values.append([rec, flow_value])
            flow_dicts[rec] = flow_dict
        except:
            print("node " + rec + " is not in the graph")

    df = pd.DataFrame(list(map(lambda val: val[1], flow_values)), columns=["DSA"])
    # df.DSA = df.DSA * scale_factor
    df.index = list(map(lambda val: val[0], flow_values))
    df["Recp"] = df.index
    gd = graphs_dict(flow_dicts, gpf)

    return df, gd




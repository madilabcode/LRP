from cmath import isnan
from networkx.testing.test import run
import pandas as pd
import numpy as np
import math
from rpy2.robjects import NULL, r
# import matplotlib.pyplot as plt
from sklearn.metrics import roc_curve, roc_auc_score
from sklearn.ensemble import RandomForestClassifier
import scipy.stats as sc
import networkx as nx
import pickle
from Codes import CERNO as ce
from Codes import utils as utils
from sklearn.metrics import mutual_info_score as mis
from Codes import pyDictToR as pyd
import concurrent.futures
import warnings
from scipy.stats import ranksums
import statsmodels.stats.multitest
import dill
warnings.filterwarnings('ignore') 


CUTOFF = 0.0
PA = pd.read_csv("./files/humanProtinAction.csv")
class Obj_dict:
    def __init__(self):
        self.dic = {}

    def add_obj(self, key, obj):
        self.dic[key] = obj

    def return_GD_from_key(self, key):
        return self.dic[key]

    def save_obj(self, name):
        utils.save_obj(self, name)


class graphs_dict:
    def __init__(self, flow_dics, capcity_network,do_permutation = True):
        self.flow_dics = flow_dics
        self.capcity_network = capcity_network
        self.max_flow_dict = self.calculate_max_flow_for_all_recp()
        self.max_multy_flow = None
        if do_permutation:
             self.perm_p_values =  self.perm_test_for_all_recptors()

    def sub_capcity_gaph_to_dict(self):
        sub = nx.to_pandas_edgelist(self.capcity_network)
        sub = sub.loc[sub["capacity"] < np.inf, :]
        gp = nx.from_pandas_edgelist(sub, "Source", "Target", edge_attr="capacity", create_using=nx.DiGraph())
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
        return graphs_dict(self.flow_dics, self.capcity_network,do_permutation=False)

    def calculate_max_flow(self, recp):
        dfFlow = self.make_df_flow(recp)
        max_flow = 0

        for row in dfFlow.iterrows():
            row = row[1]
            if row["Target"] == "sink" :
                max_flow += row["flow"]

        return max_flow

    def calculate_max_flow_for_all_recp(self):
        max_flow_dict = {}
        for recp in self.flow_dics.keys():
            max_flow_dict[recp] = self.calculate_max_flow(recp)

        return max_flow_dict
    
    def single_recptor_flow(self,recp):
        graph = self.make_df_flow(recp=recp)
        graph = graph.loc[graph.Target !="sink"]
        return graph 
    
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
        df.columns = ["Source", "Target", "flow"]
        return df

    def Tfs_Per_Recp_in_max_flow(self):
        tf_rec = {}
        for rec in self.flow_dics.keys():
            df = self.make_df_flow(rec)
            df = df.loc[df["Target"] == "sink", :]
            tf_rec[rec] = len(df["Source"])

        return tf_rec

    @staticmethod
    def perm_for_tf(cap):
        cap_df = graphs_dict.premutate_sub_graph(cap)
        graph =  nx.from_pandas_edgelist(cap_df, "Source", "Target", edge_attr="capacity", create_using=nx.DiGraph())
        flow_value,_ =  nx.maximum_flow(graph, "source_node" ,"sink", flow_func=nx.algorithms.flow.dinitz)
        return flow_value
   
    @staticmethod
    def flow_to_single_tf(args,num_of_perm=10):
        pa,tf,remove_tfs = args
        remove_tfs.remove(tf)
        pa = pa.loc[~(pa.Target.isin(remove_tfs)) & ~(pa.Source.isin(remove_tfs))]
        print(tf)
        gpf = nx.from_pandas_edgelist(pa, "Source", "Target", edge_attr="capacity", create_using=nx.DiGraph())
        orig_flow,_ = nx.maximum_flow(gpf, "source_node", tf, flow_func=nx.algorithms.flow.dinitz)
        flows =  np.array([graphs_dict.perm_for_tf(pa.copy()) for i in range(num_of_perm)])
        p_value = (flows >= orig_flow).sum() / len(flows)  
        return tf, p_value
    
    def flow_dict_to_single_tf(self,tf,pa=None):
        if pa is None:
            pa = self.build_multy_source_network()
            pa.columns = ["Source","Target","capacity"]
        remove_tfs = list(np.unique(list(pa.loc[pa["Target"] == "sink", "Source"])))        
        remove_tfs.remove(tf)
        pa = pa.loc[~(pa.Target.isin(remove_tfs)) & ~(pa.Source.isin(remove_tfs))]
        gpf = nx.from_pandas_edgelist(pa, "Source", "Target", edge_attr="capacity", create_using=nx.DiGraph())
        _,flow_dict = nx.maximum_flow(gpf, "source_node", tf, flow_func=nx.algorithms.flow.dinitz)
        df = self.make_df_flow(flow=flow_dict)
        df = df.loc[(df.Source!="source_node") & (df.Target !="sink")]
        return  df
    
    def flow_to_tfs(self,pa=None):
        if pa is None:
            pa = self.run_multiy_source_flow()
        all_tfs = np.unique(list(pa.loc[pa["Target"] == "sink", "Source"]))
        tf_flow = {}
        with concurrent.futures.ProcessPoolExecutor(max_workers=3) as executor:
            results = list(executor.map(self.flow_to_single_tf,
                                        [(pa.copy(), tf, list(all_tfs).copy()) for tf in all_tfs]))
        for element in results:
            tf_flow[element[0]] = element[1]
           
        tf_effects = pd.DataFrame({"flow":tf_flow.values()}, index =tf_flow.keys())
        tf_effects = tf_effects.sort_values(by="flow", axis=0, ascending=False)
        return tf_effects
    
    def flow_to_all_tf(self,pa=None):
        if pa is None:
            pa = self.build_multy_source_network()
        pa.columns = ["Source","Target","capacity"]
        gpf = nx.from_pandas_edgelist(pa, "Source", "Target", edge_attr="capacity", create_using=nx.DiGraph())
        cent = nx.degree_centrality(gpf)
        max_flow, flow_network = nx.maximum_flow(gpf, "source_node","sink", flow_func=nx.algorithms.flow.dinitz)
        flow_network = self.make_df_flow(flow=flow_network)
        
        all_tfs = flow_network.loc[flow_network["Target"] == "sink", ["Source","flow"]]
        all_tfs = all_tfs.groupby("Source")["flow"].sum()
        return  pa, pd.DataFrame(all_tfs.sort_values(axis=0, ascending=False))


    def calculate_significant_tf(self, recp):
        df = self.make_df_flow(recp)
        all_tf = np.unique(df.loc[df["Target"] == "sink", "Source"])
        tf_effects = {"tfs": [], "flpw_delta": []}

        for tf in all_tf:
            pa = nx.to_pandas_edgelist(self.capcity_network)
            pa.columns = ["Source", "Target", "capacity"]
            pa["isNotTf"] = pa.apply(
                lambda row: True if str(row["Source"]).find(tf) < 0 and str(row["Target"]).find(tf) < 0 else False,
                axis=1)
            pa = pa.loc[pa["isNotTf"], :]
            gpf = nx.from_pandas_edgelist(pa, "Source", "Target", edge_attr="capacity", create_using=nx.DiGraph())
            flow_value = nx.maximum_flow(gpf, recp, "sink", flow_func=nx.algorithms.flow.dinitz)[0]
            tf_effects["tfs"].append(tf)
            tf_effects["flow_delta"].append(self.max_flow_dict[recp] - flow_value)

        tf_effects = pd.DataFrame.from_dict(tf_effects)
        tf_effects = tf_effects.sort_values(by="flow_delta", axis=0, ascending=False)
        return tf_effects

    @staticmethod
    def _calculate_flow_delta_for_df(args):
        gene, pat, max_flow = args
        if gene != "source_node" and gene != "sink":
            pat["isNotTf"] = pat.apply(
                lambda row: True if row["source"] != gene  and row["target"] != gene else False, axis=1)
            pat = pat.loc[pat["isNotTf"], :]
            gpf = nx.from_pandas_edgelist(pat, "Source", "Target", edge_attr="capacity", create_using=nx.DiGraph())
            flow_value = nx.maximum_flow(gpf, "source_node", "sink", flow_func=nx.algorithms.flow.dinitz)[0]
            return gene, max_flow - flow_value

    def build_multy_source_network(self):
        pa = nx.to_pandas_edgelist(self.capcity_network)
        pa.columns = ["Source", "Target", "capacity"]
    
        for recp in self.flow_dics.keys():
            pa.loc[pa.shape[0]] = ["source_node",recp,math.inf]
        
        return pa 

    def run_multiy_source_flow(self):
        pa = self.build_multy_source_network()
        gpf = nx.from_pandas_edgelist(pa, "Source", "Target", edge_attr="capacity", create_using=nx.DiGraph())
        max_flow = nx.maximum_flow(gpf, "source_node","sink", flow_func=nx.algorithms.flow.dinitz)[0]
        self.max_multy_flow = max_flow
        return pa 

    def calculate_significant_tf_Multy_sinc(self):
        pa = self.run_multiy_source_flow()
        all_genes = np.unique(list(pa.loc[pa["Target"] == "sink", "Source"]))
        delta_dict = {}
        with concurrent.futures.ProcessPoolExecutor(max_workers=3) as executor:
            results = list(executor.map(self._calculate_flow_delta_for_df,
                                        [(gene, pa.copy(), self.max_multy_flow) for gene in all_genes]))
        for element in results:
            delta_dict[element[0]] = element[1]
           
        tf_effects = pd.DataFrame({"flow":delta_dict.values()}, index =delta_dict.keys())
        tf_effects = tf_effects.sort_values(by="flow", axis=0, ascending=False)
        return tf_effects


    @staticmethod
    def premutate_sub_graph(sub_graph,t=10):
        inf_sub_graph = sub_graph.loc[sub_graph.capacity == np.inf,:]
        sub_graph = sub_graph.loc[sub_graph.capacity < np.inf,:]

        for _ in range(t):
            sub_graph["Target"] = np.random.permutation(sub_graph["Target"])
        self_circ_rows = sub_graph.loc[sub_graph["Target"] == sub_graph["Source"]].index

        for row in self_circ_rows:
            target = sub_graph.loc[row]["Target"]
            perm_row = sub_graph.loc[(sub_graph["Target"] != target) & (sub_graph["Source"] != target) ,:].sample()
            sub_graph.loc[row,"Target"] = perm_row["Target"].to_list()[0]
            sub_graph.loc[perm_row.index,"Target"] = target
        return pd.concat([sub_graph,inf_sub_graph])
        
    def edge_dgree_perm(self, source_node,cap=None):
        if cap is None:
            cap = self.capcity_network
        
        cap_df = nx.to_pandas_edgelist(cap)
        cap_df.columns = ["Source", "Target", "capacity"]
        ##### 3 bins
        cap_df.sort_values("capacity",inplace=True)
        top1 = cap_df.iloc[int(cap_df.shape[0]/3)]["capacity"]
        top2 = cap_df.iloc[int(2*cap_df.shape[0]/3)]["capacity"]
        sub_cup1 = self.premutate_sub_graph(cap_df.loc[cap_df.capacity < top1,:])
        sub_cup2 = self.premutate_sub_graph(cap_df.loc[(cap_df.capacity >= top1) & (cap_df.capacity <= top2) ,:])
        sub_cup3 = self.premutate_sub_graph(cap_df.loc[cap_df.capacity > top2,:])
        cap_df = pd.concat([sub_cup1,sub_cup2,sub_cup3])
        ###on bine!
        cap_df = self.premutate_sub_graph(cap_df)
        graph =  nx.from_pandas_edgelist(cap_df, "Source", "Target", edge_attr="capacity", create_using=nx.DiGraph())
        flow_value, flow =  nx.maximum_flow(graph, source_node ,"sink", flow_func=nx.algorithms.flow.dinitz)
        return flow_value

    def calculate_p_value_for_recp(self, rec, num_of_perm=100):
        print(rec)
        orig_flow, _ = nx.maximum_flow(self.capcity_network, rec ,"sink", flow_func=nx.algorithms.flow.dinitz)
        flows =  np.array([self.edge_dgree_perm(rec) for i in range(num_of_perm)])
        p_value = (flows >= orig_flow).sum() / len(flows)  
        print(p_value)
        return p_value
    
    def calculate_p_value_for_tfs(self,num_of_perm=100):
        pa, flows = self.flow_to_tfs() 
        pa.columns = ["Source","Target","capacity"]
    
        for i in range(num_of_perm):
            perm_pa = self.premutate_sub_graph(pa)
            if i == 0:
                df = self.flow_to_tfs(perm_pa)[1]
            else :
                df = df.merge(self.flow_to_all_tf(perm_pa)[1],left_index=True, right_index=True)
        
        flows = pd.DataFrame(flows)
        flows["gene"] = flows.index
        flows = flows.loc[flows.index.isin(df.index)]
        flows["p_value"] = flows.apply(lambda x: (df.loc[x[1]] > x[0]).sum() / num_of_perm ,axis=1)
        return  flows[["flow","p_value"]]
            
    def perm_test_for_all_recptors(self):
        pvals = {}
        for recp in self.flow_dics.keys():
            pvals[recp] = (self.calculate_p_value_for_recp(recp))
        return pvals
    
    def max_flow_withno_recptor(self,recp,exp):
        cap = self.capcity_network
        cap_df = nx.to_pandas_edgelist(cap)
        cap_df.columns = ["Source", "Target", "capacity"]
        cap_df = cap_df.loc[cap_df.Source != recp,:]
        cap_df = pd.concat([cap_df,self.pa.loc[self.pa.Source == recp,:]])
        dist_args =  utils.load_obj(r"./outputObj/dist_norm_args")
        cap_df.loc[cap_df.Source == recp,"capacity"] = 1
        graph =  nx.from_pandas_edgelist(cap_df, "Source", "Target", edge_attr="capacity", create_using=nx.DiGraph())
        flow_value, _ =  nx.maximum_flow(graph, recp ,"sink", flow_func=nx.algorithms.flow.dinitz)
        return flow_value


def format_Gene_dict(gene, from_dict=False):
    if from_dict:
        gene = gene[1:len(gene) - 1]

    fl = gene[0]
    sl = gene[1:]
    sl = str.lower(sl)

    return fl + sl

def save_obj(obj, name):
    with open(name + '.pkl', 'wb') as f:
        dill.dump(obj, f, dill.HIGHEST_PROTOCOL)

def load_obj(name):
    with open(name + '.pkl', 'rb') as f:
        return dill.load(f)

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

def mi(exp,pa, recps_for_roc=None):
    clusterex = exp.copy()
    clusterex["not_zero"] = clusterex.apply(lambda x: x.astype(bool).sum(), axis=1)
    if recps_for_roc is None:
         clusterex = clusterex.loc[clusterex.not_zero > clusterex.shape[1] * CUTOFF, :]
    else:
        clusterex = clusterex.loc[(clusterex.not_zero > clusterex.shape[1] * CUTOFF) | (clusterex.index.isin(recps_for_roc)) , :]

    clusterex.drop("not_zero",axis=1,inplace = True)

    pa = pa.loc[(pa.Source.isin(exp.index)) &  (pa.Target.isin(exp.index) )]
    pa["wights"] = pa.apply(lambda x : mis(exp.loc[x.Source],exp.loc[x.Target]),axis=1)
    pa["wights"] = pa.apply(lambda x: x.wights if x.Source in clusterex.index and x.Target in clusterex.index else 0 ,axis=1)
    pa["wights"] = (pa["wights"] - pa["wights"].min())/(pa["wights"].max() - pa["wights"].min())
    pa["capacity"] = pa["wights"]
    pa.drop("wights",axis=1,inplace = True)

    return pa

def build_flowing_network_with_normlized_wights(ProtAction, tfs_enrc, exp, recps_for_roc = None,wights_flag=True):
    pa = ProtAction
    pa["Source"] = pa["Input-node Gene Symbol"]
    pa["Target"] = pa["Output-node Gene Symbol"]
    pa["capacity"] = pa["Edge direction score"]
    keeps = ["Source", "Target", "capacity"]
    pa = pa[keeps]

    pa["Source"] = pa["Source"].apply(format_Gene_dict).astype(str)
    pa["Target"] = pa["Target"].apply(format_Gene_dict).astype(str)
    
    if wights_flag:
        pa = mi(exp,pa)
    pa.reset_index(drop=True,inplace=True)

    for tf in tfs_enrc:
        if tf in pa["Target"].to_list():
            pa.loc[pa.shape[0]] = [tf, "sink", np.inf]
 

    gp = nx.from_pandas_edgelist(pa, "Source", "Target", edge_attr="capacity", create_using=nx.DiGraph())
    utils.save_obj(gp, "temp_grpah")

    return gp

def hypergeomtric_test(gs, markers,genes):
    N = len(genes)
    m = len(list(filter(lambda x: x in genes,gs)))
    k = len(markers)
    x = len(list(filter(lambda x: x in markers.index,gs)))
    return 1 - sc.hypergeom.cdf(x,N,m,k) 

def wilcoxon_enrcment_test(tf,gene_list,exp):
    try:
        gene_exp = exp.loc[exp.index.isin(gene_list)]
        backround_exp = exp.loc[~exp.index.isin(gene_list)]
        return ranksums(backround_exp,gene_exp,alternative="two-sided")[1]
    except:
        print(f"problem accure in {tf}")
        return 1
                                
def enrch_tfs(exp,tfs, reduce):
    expamount = exp.apply(lambda x: x.astype(bool).sum(), axis=1)
    expamount = expamount[expamount > exp.shape[1]*0.1]
    tfs = {key: value for key,value in tfs.items() if key in expamount.index}
    tfs_scores = {tf:wilcoxon_enrcment_test(tf,tfs[tf],exp.loc[expamount.index].mean(axis=1)) for tf in tfs.keys() if len(tfs[tf])>1}
    tfs_scores = {key: value for key,value in tfs_scores.items() if not np.isnan(value)}
    values = statsmodels.stats.multitest.fdrcorrection(np.array(list(tfs_scores.values())))[1]
    tfs_scores = {key:values[index] for index,key in enumerate(tfs_scores.keys())}
    if reduce:
        tfs_scores = {tf: tfs_scores[tf] for tf in tfs_scores.keys() if tfs_scores[tf]<=0.05}
    return tfs_scores

def DSA_anaylsis(exp, recptors, ProtAction, ProtInfo, tfs,markers=None,recps_for_roc = None,reduce=True,do_permutation=True,wights_flag=True):
    tfs_scores = enrch_tfs(exp,tfs,reduce)
    gpf = build_flowing_network_with_normlized_wights(ProtAction, tfs_scores, exp,recps_for_roc=recps_for_roc,wights_flag=wights_flag)
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
    gd = graphs_dict(flow_dicts, gpf,do_permutation)
    return df, gd

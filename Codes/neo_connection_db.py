from py2neo import Graph, Node, Relationship, NodeMatcher, RelationshipMatcher, Neo4jError
import numpy as np
import pandas as pd
import json
import socket
from Codes import Tf_Graph as tf
from Codes import utils


def table_to_dict(df, key1, key2):
    keys = np.unique(df[key1])
    dic = {}
    for key in keys:
        values = list(np.unique(df.loc[df[key1] == key, key2]))
        dic[key] = values

    return dic


def dict_to_table(dic, keysName, valuesName):
    keys = []
    vlaues = []

    for key in dic.keys():
        for vlaue in dic[key]:
            keys.append(key)
            vlaues.append(vlaue)

    df = pd.DataFrame([keys, vlaues]).transpose()
    df.columns = [keysName, valuesName]
    return df


def start_connection(userg="neo4j", passwordg="neo4j_new"):
    try:
        graph = Graph(user=userg, password=passwordg)
        return graph
    except:
        try:
            graph = Graph(user=userg, password=passwordg)
            return graph
        except socket.error as e:
            print("Could not oonnect to graph - return error " + str(e))
            return None


def ceack_if_node_in_graph(graph, node):
    props = []
    keys = []

    for prop in node.keys():
        props.append(prop)
        keys.append(node[prop])

    return return_node_by_list_of_prop(graph, props, keys)


def add_node_to_neo(graph, ClusterName, treatment, proj="Michael_proj"):
    try:
        type_of_cell, subTypeCell = return_type_and_subtype(ClusterName)
    except:
        return Exception("no data found")
    node = Node("Cell", Type=type_of_cell, SubType=subTypeCell, Cluster=ClusterName,
                treatment=treatment, project=proj)

    exnode = ceack_if_node_in_graph(graph, node)
    if len(exnode) == 0:
        graph.create(node)
        nodes = return_list_of_nodes_by_prop(graph, "type", type_of_cell)
        for tnode in nodes:
            if tnode != node:
                make_smae_type_relationship(graph, node, tnode, type_of_cell)
        return node
    else:
        return exnode


def add_signal_relationship_to_neo(graph, sourceNode, TargetNode, DSA_mean, DSA_table, legRet, graphObj,
                                   DE_interaction=None):
    legRetDict = table_to_dict(legRet, "Ligand", "Receptor")
    graph.create(Relationship(sourceNode, "Signal", TargetNode, DSA_mean=DSA_mean, Recp=list(DSA_table["Recp"]),
                              DSA_Score=list(DSA_table["DSA"]) \
                              , DE_Recp=None, LegRec=json.dumps(legRetDict),
                              capcity_network=utils.json_zip(json.dumps(graphObj.sub_capcity_gaph_to_dict())),
                              DE_LegRec=None))


def crate_all_from_nodes(graph, toNode, DSA_means, DSA_Tabls, legRets, ObjDictFrom, treatment, proj="Michael_proj"):
    DSA_Tabls = dict(DSA_Tabls)
    legRets = dict(legRets)
    DSA_means = table_to_dict(DSA_means, "Cluster", "DSA_Mean")

    nodes = {}

    for cluster in DSA_Tabls.keys():

        try:
            nodes[cluster] = \
                return_node_by_list_of_prop(graph, ['Cluster', 'treatment', 'project'], [cluster, treatment, proj])[0]
        except:
            try:
                nodes[cluster] = add_node_to_neo(graph, cluster, treatment, proj)
            except:
                print(cluster + " does not have all requride info")

    for cluster in nodes.keys():
        try:
            add_signal_relationship_to_neo(graph, nodes[cluster], toNode, float(DSA_means[cluster][0]),
                                           DSA_Tabls[cluster],
                                           legRets[cluster], ObjDictFrom.return_GD_from_key(cluster).update_obj())
        except:
            print(cluster + " does not have all requride info2")
    return nodes


def delete_nodes_by_prop(graph, prop=None, key=None, delete_all_grpah=True):
    matcher = NodeMatcher(graph)
    node_it = matcher.match("Cell").all()
    subNode = [node for node in node_it if node[prop] == key or delete_all_grpah]

    matcher = RelationshipMatcher(graph)
    rel_int = matcher.match().all()
    rel_int = list(
        filter(lambda rel: rel.start_node[prop] == key or rel.end_node[prop] == key or delete_all_grpah, rel_int))

    for node in subNode:
        graph.delete(node)

    for rel in rel_int:
        graph.delete(rel)


def return_list_of_nodes_by_prop(graph, prop, key):
    matcher = NodeMatcher(graph)
    node_it = matcher.match("Cell").all()
    nodes = [node for node in node_it if node[prop] == key]
    return nodes


def compere_props_of_entty(ent, props, keys):
    for index, prop in enumerate(props):
        if not ent[prop] == keys[index]:
            return False
    return True


def return_node_by_list_of_prop(graph, props, keys):
    matcher = NodeMatcher(graph)
    node_it = matcher.match("Cell").all()
    nodes = [node for node in node_it if compere_props_of_entty(node, props, keys)]
    return nodes


def make_cluster_trat_relationship(graph, node1, node2, DE_interaction):
    De_dict = {key: (list(value["feture"]), list(value["importances value"])) for key, value in DE_interaction.items()}
    graph.create(
            Relationship(node1, "Same_Cluster", node2, DE_LegRec=json.dumps(De_dict)))


def make_smae_type_relationship(graph, node1, node2, CellType):
    graph.create(Relationship(node1, "same_cell_type", node2, CellsType=CellType))
    graph.create(Relationship(node2, "same_cell_type", node1, CellsType=CellType))


def find_relationship_between_nodes(graph, node1, node2):
    matcher = RelationshipMatcher(graph)
    return list(matcher.match([node1, node2]))


def find_relationships_by_type_and_props(graph, typeRe, props=None, keys=None):
    matcher = RelationshipMatcher(graph)
    rel_int = matcher.match(r_type=typeRe).all()

    if props is not None:
        rels = [rel for rel in rel_int if compere_props_of_entty(rel, props, keys)]
        return rels
    return list(rel_int)


def make_type_and_subtype(clusters, proj):
    types = {}
    subtypes = {}

    for cluster in clusters:
        types[cluster] = input("Enter the Cell Type of cluster " + str(cluster) + " ? ")
        subtypes[cluster] = input(
            "Enter the Cell sub Type of cluster " + str(cluster) + " ?, if there is not enter None ")

    print(proj)
    tf.save_obj(types, f"outputObj/{proj}_types")
    tf.save_obj(subtypes, f"outputObj/{proj}_subtypes")

    return types, subtypes


def return_type_and_subtype(cluster):
    try:
        annot = pd.read_csv("Annot.csv")
    except:
        raise Exception("No anontaion for type and sub type found")

    typeCell = list(annot.loc[annot["Cluster"].astype(str) == str(cluster), "Type"])[0]
    subTypeCell = list(annot.loc[annot["Cluster"].astype(str) == str(cluster), "subType"])[0]

    return typeCell, subTypeCell


def return_all_DSA_Mean_of_proj(graph, proj, Treat=None):
    matcher = RelationshipMatcher(graph)
    rel_int = list(matcher.match(r_type="Signal").all())
    if Treat is None:
        rel_int = list(
            filter(lambda rel: rel.end_node["project"] == proj and rel.start_node["project"] == proj, rel_int))
    else:
        rel_int = list(filter(lambda rel: rel.end_node["project"] == proj and rel.start_node["project"] == proj and
                                          rel.end_node["treatment"] == Treat and rel.start_node["treatment"] == Treat,
                              rel_int))

    sourceNodes = [f"{rel.start_node['Type']}/{rel.start_node['SubType']}/{rel.start_node['Cluster']}"
                   if rel.start_node['SubType'] != "None" else f"{rel.start_node['Type']}/{rel.start_node['Cluster']}"
    if rel.start_node['Type'] != rel.start_node['Cluster'] else rel.start_node['Cluster'] for rel in rel_int]

    targetNodes = [f"{rel.end_node['Type']}/{rel.end_node['SubType']}/{rel.end_node['Cluster']}"
                   if rel.end_node['SubType'] != "None" else f"{rel.end_node['Type']}/{rel.end_node['Cluster']}"
    if rel.end_node['Type'] != rel.end_node['Cluster'] else rel.end_node['Cluster'] for rel in rel_int]

    Scores = [rel["DSA_mean"] for rel in rel_int]
    df = pd.DataFrame([sourceNodes, targetNodes, Scores]).transpose()
    df.columns = ["From Cluster", "To Cluster", "DSA-MEAN"]
    return df


def return_legRet_and_DSA_of_proj_and_Treat(graph, proj, Treat):
    matcher = RelationshipMatcher(graph)
    rel_int = list(matcher.match(r_type="Signal").all())
    DSA_Tables = {}
    legRet = {}

    rel_int = list(filter(lambda rel: rel.end_node["project"] == proj and rel.start_node["project"] == proj and
                                      rel.end_node["treatment"] == Treat and rel.start_node["treatment"] == Treat,
                          rel_int))

    clusters = np.unique([rel.end_node["Cluster"] for rel in rel_int])

    for fromCluster in clusters:
        DSA_Tables[fromCluster] = {}
        legRet[fromCluster] = {}

        for toCluster in clusters:
            if fromCluster != toCluster:
                sig = list(filter(
                    lambda rel: rel.start_node["Cluster"] == fromCluster and rel.end_node["Cluster"] == toCluster,
                    rel_int))
                if len(sig) < 1:
                    print(
                        f"There is no edge for signal: {fromCluster} to {toCluster} for Project {proj} and Treatment {Treat}")
                    continue
                elif len(sig) > 1:
                    raise Exception(
                        f"There is more then one edge for signal: {fromCluster} to {toCluster} for Project {proj} and Treatment {Treat}")

                DSA_Tables[fromCluster][toCluster] = [[recp for recp in sig[0]["Recp"]],
                                                      [score for score in sig[0]["DSA_Score"]]]
                legRet[fromCluster][toCluster] = dict_to_table(json.loads(sig[0]["LegRec"]), "Ligand", "Receptor")

    for fromCluster in DSA_Tables.keys():
        for toCluster in DSA_Tables[fromCluster].keys():
            if fromCluster != toCluster:
                DSA_Tables[fromCluster][toCluster] = pd.DataFrame(DSA_Tables[fromCluster][toCluster]).transpose()
                DSA_Tables[fromCluster][toCluster].columns = ["Recp", "DSA"]

    return DSA_Tables, legRet


def return_all_DSA_and_legRet_of_proj(graph, proj, Treats):
    dic = {treat: return_legRet_and_DSA_of_proj_and_Treat(graph, proj, treat) for treat in Treats}
    DSA_Tables = {treat: dic[treat][0] for treat in dic.keys()}
    legRet = {treat: dic[treat][1] for treat in dic.keys()}
    return DSA_Tables, legRet


def return_rel_capcity_network(graph, proj, Treat, toCluster, fromCluster):
    toNode = return_node_by_list_of_prop(graph, ["Cluster", "treatment", "project"], [toCluster, Treat, proj])[0]
    fromNode = return_node_by_list_of_prop(graph, ["Cluster", "treatment", "project"], [fromCluster, Treat, proj])[0]
    sig = list(find_relationship_between_nodes(graph, fromNode, toNode))
    sig = list(filter(lambda s: str(s).find("Signal") >= 0, sig))[0]
    return tf.Obj_dict.make_obj_out_of_capcity_network(json.loads(utils.json_unzip(sig["capcity_network"])),
                                                       sig["Recp"])

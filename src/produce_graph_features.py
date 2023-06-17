import os
import re
import glob 
import argparse
from typing import Dict
from pathlib import Path

import numpy as np
import pandas as pd

import igraph as ig

# # from code.parser_utils import parse_all
# from graph_features.code.parser_utils import parse_all

##############################################################
##############################################################
##############################################################
# parse output of rinerator

def parse_grafo(prot="7KVE_relaxed", path_data="../data/igraph/"):
    
    grafo_path = os.path.join(path_data, "in", f'{prot}_h.sif')
    df_grafo = pd.read_csv(grafo_path, sep="\s+|:", header=None, engine="python")

    dic_rename = {1: "Pos1", 3: "Res1", 
                  4: "inter1", 5: "inter2",
                  7: "Pos2", 9: "Res2"}

    df_grafo = df_grafo.rename(columns=dic_rename)

    df_grafo = df_grafo[list(dic_rename.values())]

    return df_grafo


def parse_nint(prot="7KVE_relaxed", path_data="../data/igraph/"):

    nint_file = os.path.join(path_data, "in", f"{prot}_h_nrint.ea")
    
    df_nint = pd.read_csv(nint_file, sep="\s+|:", skiprows=1, header=None, engine="python")

    dic_rename = {1: "Pos1", 3: "Res1", 
                  4: "inter1", 5: "inter2",
                  7: "Pos2", 9: "Res2", 
                  10: "Nint"}

    df_nint = df_nint.rename(columns=dic_rename)

    df_nint = df_nint[list(dic_rename.values())]

    df_nint["inter1"] = df_nint["inter1"].apply(lambda x: x.strip("("))
    df_nint["inter2"] = df_nint["inter2"].apply(lambda x: x.strip(")"))

    return df_nint


def parse_scores(prot="7KVE_relaxed", path_data="../data/igraph/"):
    
    scores_path = os.path.join(path_data, "in", f"{prot}_h_intsc.ea")
    df_scores = pd.read_csv(scores_path, sep="\s+|:", skiprows=1, header=None, engine="python")

    dic_rename = {1: "Pos1", 3: "Res1", 
                  4: "inter1", 5: "inter2",
                  7: "Pos2", 9: "Res2",
                  10: "Scores"}

    df_scores = df_scores.rename(columns=dic_rename)

    df_scores = df_scores[list(dic_rename.values())]

    df_scores["inter1"] = df_scores["inter1"].apply(lambda x: x.strip("("))
    df_scores["inter2"] = df_scores["inter2"].apply(lambda x: x.strip(")"))

    return df_scores

def concat_to_order(row):
    return "|".join(sorted([str(row[col]) for col in "Pos1 Pos2 Res1 Res2".split()]))

def parse_all(prot="7KVE_relaxed", path_data="../data/igraph/", unidirectional=True, filter_inter=False):
    
    df_grafo = parse_grafo(prot, path_data)
    df_nint = parse_nint(prot, path_data)
    df_scores = parse_scores(prot, path_data)

    if (df_grafo.shape[0] == df_nint.shape[0]) and (df_grafo.shape[0] == df_scores.shape[0]):

        keys = df_grafo.columns.tolist()

        df = df_grafo.merge(df_nint, on=keys, how="outer").merge(df_scores, on=keys, how="outer")
        
        if df_grafo.shape[0] == df.shape[0]:
            
            if filter_inter:
                
                inter1_list = df["inter1"].unique()

                print("\nAvailable interactions (type 1):")

                for i, x in enumerate(inter1_list):
                    print(f'{i+1} - {x}')

                inter1 = int(input("\nPlease type the corresponding number of the desired interaction in the list above: "))
                inter1 = inter1_list[inter1-1]

                #################################################
                print("\n")
                print("="*50)
                print("\n")
                #################################################

                inter2_list = df[df["inter1"] == inter1]["inter2"].unique()

                print("Available interactions (type 2):")

                for i, x in enumerate(inter2_list):
                    print(f'{i+1} - {x}')

                inter2 = int(input("\nPlease type the corresponding number of the desired interaction in the list above: "))
                inter2 = inter2_list[inter2-1]
                
                df = df[(df["inter1"] == inter1) & (df["inter2"] == inter2)].copy()
                
                str_inter = f"_{inter1}_{inter2}"
                
            else:
                str_inter = ""

            ################################
            
            # maybe this is unnecessary in the cases that the filter above is applied
            # but let's do it either way -- wont harm.
            
            if unidirectional:
                
                # the graph must be unidirectional. 
                # in order to eliminate possible bi-directional edges, we create an 
                # auxiliary column with a sorted string of the first columns.
                # example, for the data below:

                # Pos1	 Pos2	Res1	Res2	Nint	Scores
                # 305	 306	THR   	ARG 	1   	0.0027
                # 306	 305	ARG   	THR 	2   	0.1327

                # we have:

                # Pos1	 Pos2	Res1	Res2	Nint	Scores    order
                # 305	 306	THR   	ARG 	1   	0.0027    305|306|ARG|THR
                # 306	 305	ARG   	THR 	2   	0.1327    305|306|ARG|THR

                # by using this new column, we can detect that both rows describe an edge between the same 
                # nodes, thus, only one of them must be chosen (the one with smallest score, see below)

                df["order"] = df.apply(concat_to_order, axis=1)

                # we keep only the smallest score
                # (this idea came up in the call with tiago on 14/aug/2021: possible new contribution)
                # notice: i get the index of the minimum score among each group, and then
                # use this index to get only this single edge within each group ;)

                df = df.loc[df.groupby("order")["Scores"].idxmin()].copy()
                
                # to identify the final file
                str_uni = "unidirectional"
                
            else:
                
                # here i only consider the smallest score, but keep any possible bi-directed edge
                
                df = df.groupby("Pos1 Pos2 Res1 Res2".split()).min().reset_index(drop=False).copy()
            
                str_uni = "bidirectional"

            ######################################
            
            print("\nDistribution of interaction (type 1):")
            #display(df["inter1"].value_counts())
            print(df["inter1"].value_counts())
            
            print("\nDistribution of interaction (type 2):")
            #display(df["inter2"].value_counts())
            print(df["inter2"].value_counts())
            
            # selecting cols of the final dataframe
            # df = df["Pos1 Pos2 Res1 Res2 Nint Scores".split()].copy()
            # new: let's export all columns, to put in the supplementary tables, that's why it's commented above
            
#             data_path = os.path.join(path_data, 'results', prot)
#             Path(data_path).mkdir(parents=True, exist_ok=True) 

            path_graph_data = os.path.join(path_data, "out", f"{prot}_graph_data_{str_uni}{str_inter}.csv")
            df.to_csv(path_graph_data, index=False)
            print(f"\nParser done and graph data exported to {path_graph_data}!")
            
        else:
            print("ERROR! REVIEW FILES/CODE")
    else:
        print("ERROR! REVIEW FILES/CODE")

    return df, str_uni, str_inter

##############################################################
##############################################################
##############################################################
# produce graph features

def make_simple(g):
    '''
    checks if graph is simples. If not, simplify it, by combining 
    egdes using the minimum value
    '''

    if g.is_simple():

        print("\nGraph is simple (no loops or multiple edges)")

    else:

        print("\nGraph was not originally simple (no loops or multiple edges)")

        g = g.simplify(combine_edges="min")

        print("Graph was simplified!")
        
    return g


#_________________________________________________


def check_number_of_features(degree, betweenness, closeness, burt, authority, page_rank, k_core):
    '''
    checks if the number of calculated graph features is consistent for all features
    if not, raises an error
    '''

    check_len = [degree, betweenness, closeness, burt, authority, page_rank, k_core]

    n = len(check_len[0])

    if not all(len(x) == n for x in check_len):

        raise ValueError("There's a problem with the number of features!\nCheck code!")
        
    else:
        
        print("\nEverythings all right with the number of features!")


#_________________________________________________


def extract_graph_features(prot="7KVE_relaxed", path_data="../data/igraph/",
                           str_uni="unidirectional", unidirectional=True,
                           str_inter="", filter_inter=False, 
                           write=True,
                           weighted=False):
    '''
    this functions calculates the graph features, readin the .csv file generated by the
    parser from the RINerator files.
    if the csv file doesn't exist yet, it calls the parser.
    
    args:
        - prot (str): string indicating the protein to be used (its code, as used in the directories and files names);
        - str_uni (str): string which indicates if the graph is uni/bidirectional. Used to identify the respective csv file;
        - unidirectional (bool): boolean flag used in the call of the parser, if csv file doesn't exist yet;
        - str_inter (str): string which indicates the kind of interaction. Used to identify the respective csv file;
        - filter_inter (bool): boolean flag used in the call of the parser, if csv file doesn't exist yet;
        - write (bool): boolean flag to control wheter or not the final dataframe of graph features is to be saved.
        # new 12/04/23 - during paper review
        - weighted (bool): if True, use "Score" as the edges weights.
    '''
    
    # if the csv files already exists
    try:
        
        file = os.path.join(path_data, 'out', f"{prot}_graph_data_{str_uni}{str_inter}.csv")
        df_parser = pd.read_csv(file)
    
    # if not, it's because the parser wasn't executed yet. So we execute it!
    except:
        
        print("\nParser will be executed now!\n")
        
        df_parser, str_uni, str_inter = parse_all(prot=prot, 
                                                  unidirectional=unidirectional, filter_inter=filter_inter,
                                                  path_data=path_data)
        
        print("\nParser executed successfully!\n")
        
    # only use these columns
    df_parser = df_parser["Pos1 Pos2 Res1 Res2 Nint Scores".split()].copy()

    if weighted:
        
        # new 12/04/23: the weights must be all positive, but we have negative interaction scores...
        # so, let's normalize them between 1 and 10
        scores = df_parser["Scores"].values
        df_parser["weight"] = ((scores - scores.min())/(scores.max() - scores.min())) * (10 - 1) + 1
        
#         display(df_parser["weight"].describe())
        
        # https://stackoverflow.com/questions/44400345/create-igraph-graph-from-pandas-dataframe
        g = ig.Graph.TupleList(df_parser["Pos1 Pos2 weight".split()].itertuples(index=False), 
                               directed=False,
#                                weights=True,
                               edge_attrs=["weight"],
                              )
        weighted_str = "_weighted"
    
    else:
        
        # it's important to model as a undirected, unweighted graph!
        # (meeting with tiago on 20/11/2021)
        g = ig.Graph.TupleList(df_parser.itertuples(index=False), directed=False, weights=False)
        weighted_str = ""
    
    g = make_simple(g)

    ###############################################

    degree = g.degree()
    
    if weighted:
        burt = g.constraint(weights="weight")
        closeness = g.closeness(weights="weight")
        betweenness = g.betweenness(weights="weight")
        authority = g.authority_score(weights="weight")
        page_rank = g.personalized_pagerank(weights="weight")
        k_core = g.coreness()
    else:
        burt = g.constraint()
        closeness = g.closeness()
        betweenness = g.betweenness()
        authority = g.authority_score()
        page_rank = g.personalized_pagerank()
        k_core = g.coreness()
    
    check_number_of_features(degree, betweenness, closeness, burt, authority, page_rank, k_core)

    ###############################################
    
    vertices = [x["name"] for x in g.vs]

    data = {"node" : vertices,
            "degree" : degree,
            "betweenness" : betweenness,
            "closeness" : closeness,
            "burt" : burt,
            "authority" : authority,
            "page_rank" : page_rank,
            "k_core" : k_core}

    df_graph = pd.DataFrame(data)

    df_graph["node"] = df_graph["node"].astype(int)

    if write:
        save_path = os.path.join(path_data, 'out', f"{prot}_graph_features_{str_uni}{str_inter}.csv")
        df_graph.to_csv(save_path, index=False)
        # print(f"\nFile saved in {save_path}!\n")
        
    return df_graph

#______________________________________________


def get_user_args():

    args = argparse.ArgumentParser()

    args.add_argument('--prot',
                      default="7KVE_relaxed", 
                      help='Protein PDB code.')

    # Add the arguments
    args.add_argument('--ip', '--input_path',
                      type=str,
                      default="../data/igraph/",
                      help='Path to the folder with contains the txt files to concat')
    
    args.add_argument('--unidirectional', 
                      default=True, 
                      help='Select only the lowest connection score.')
    args.add_argument('--filter_inter', 
                      default=False, 
                      help='Select different protein interactions.')
    args.add_argument('--write', 
                      default=True, 
                      help='Write graph features.')

    return args

#_________________________________________________

# dafuq?
def main():
    parser = get_user_args()
    args = parser.parse_args()

    _ = extract_graph_features(prot=args.prot, 
                               path_data=args.ip,
                               unidirectional=args.unidirectional, 
                               filter_inter=args.filter_inter,
                               write=args.write)

#_________________________________________________
    
if __name__ == "__main__":
    main()
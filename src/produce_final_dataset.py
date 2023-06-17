from produce_graph_features import *

def read_mutations_dataset(path_mutation="../data"):
    '''
    reads and returns the mutations dataset
    '''
    
    path = f"{path_mutation}/FV_mutations.csv"

    df_mutation = pd.read_csv(path)
    
    # new 25/12/2021: now we also bring "Protein.Change", as Tiago requested
    cols = ["Legacy.Amino.Acid", "Reported.Clinical.Severity", "Merged_Clinical_Severity", "Protein.Change"]

    return df_mutation[cols]

#______________________________________________

def merge_produce_final_graph(
    prot="7KVE_relaxed", path_data = "../data/igraph/", path_mutation="../data/",
    str_uni="unidirectional", unidirectional=True,
    str_inter="", filter_inter=False, 
    write=True, merge_how="left", weighted=False
):
    '''
    read the dataset with graph features, and merge it with the mutations dataset.
    arg "merge_how" controls if NaNs are to be dropped:
    - if "left" : keep all nodes from the mutation datasets
    - if "inner" : keep only matching nodes (the same as merging with "left" followed by .dropna()
    '''
    
    df_graph = extract_graph_features(prot=prot, path_data=path_data,
                                      str_uni=str_uni, unidirectional=unidirectional, 
                                      str_inter=str_inter, filter_inter=filter_inter, 
                                      write=write, weighted=weighted)
    
    df_mutation = read_mutations_dataset(path_mutation=path_mutation)
    
    ###############################################

    parse_res_mutation_regex = r"(?P<res_orig>[A-Za-z*]+)\d+(?P<res_mutated>[A-Za-z*]+)"

    res_change_parsed = df_mutation["Protein.Change"].str.extract(parse_res_mutation_regex, expand=False)

    df_mutation = pd.concat([df_mutation, res_change_parsed], axis=1)

    ###############################################
    
    # change to "inner" if wants to drop NaNs directly...
    df_final = df_mutation.merge(df_graph, left_on="Legacy.Amino.Acid", right_on="node", how=merge_how)
    
    # align keys, just to make it clearer
    keys = ['Legacy.Amino.Acid','node']
    df_final = df_final[keys + [x for x in df_final.columns if x not in keys]].copy()
    
    return df_final

#______________________________________________

def read_all_data(prot = "7KVE_clean",
                  str_uni = "unidirectional", unidirectional = True,
                  str_inter = "", filter_inter = False,
                  weighted=False):
    

    df_mutation_graph = merge_produce_final_graph(
        prot=prot, 
        path_data = "../data/igraph/", path_mutation="../data/",
        str_uni=str_uni, unidirectional=unidirectional,
        str_inter=str_inter, filter_inter=filter_inter,
        write=True, merge_how="inner", weighted=weighted
    )

    # print(f"\n\nShape of data: {df_mutation_graph.shape}")

    ##########################################

    # if the file doesn't exist yet, run: 
    # python parse_chimera_features.py -ip=f"data/{prot}"
    # to create the file. (See "explain" notebook in chimera folder for details)

    chimera_path = os.path.join("../data/chimera", "out", f"{prot}_structure.csv")

    # to-do: CRIAR O CSV!!!!
    df_structure = pd.read_csv(chimera_path)

    ##########################################

    conservation_path = os.path.join("../data", "FV_conservation.csv")

    df_conservation = pd.read_csv(conservation_path, sep="\t")

    ##########################################

    merge_how = "inner"

    # df_mutation_graph_structure \equiv df_mgs
    df_mgs = df_mutation_graph.merge(df_structure, left_on="Legacy.Amino.Acid", right_on="number", how=merge_how)

    # print(f"\nShape of data: {df_mgs.shape}")

    ##########################################

    bring_only = [' POS', 'SCORE']

    # df_mutation_graph_structure_conservation \equiv df_mgsc
    df_mgsc = df_mgs.merge(df_conservation[bring_only], left_on="Legacy.Amino.Acid", right_on=' POS', how=merge_how)

    keys = ['Legacy.Amino.Acid', 'node', ' POS', 'number']
    df_mgsc = df_mgsc[keys + [x for x in df_mgsc.columns if x not in keys]].copy()

    # print(f"\nShape of data: {df_mgsc.shape}, {df_structure.shape}, {df_conservation.shape}")
    
    return df_mgsc, df_structure, df_conservation

#______________________________________________

def get_user_args():

    args = argparse.ArgumentParser()

    args.add_argument('--prot',
                      default="7KVE_relaxed", 
                      help='Protein PDB code.')
    args.add_argument('--str_uni', 
                      default="unidirectional", 
                      help='Select only the lowest connection score.')
    args.add_argument('--unidirectional', 
                      default=True, 
                      help='Select only the lowest connection score.')
    args.add_argument('--filter_inter', 
                      default=False, 
                      help='Select different protein interactions.')
    args.add_argument('--str_inter', 
                      default="", 
                      help='Select different protein interactions.')
    args.add_argument('--write', 
                      default=True, 
                      help='Write graph features.')

    return args

#_________________________________________________

def main():
    parser = get_user_args()
    args = parser.parse_args()
    
    _ = read_data(
            prot=args.prot, 
            unidirectional=args.unidirectional, 
            filter_inter=args.filter_inter,
            str_uni = args.str_uni,
            str_inter = args.str_inter,
            weighted=False
    )

#_________________________________________________
    
if __name__ == "__main__":
    main()
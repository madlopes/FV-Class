import os

import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
from matplotlib.pyplot import cm
from matplotlib.colors import to_rgba
from matplotlib.patches import Rectangle
from matplotlib.collections import PathCollection
from matplotlib.legend_handler import HandlerPathCollection, HandlerLine2D

import seaborn as sns
sns.set_theme(style="darkgrid")

import ast

from sklearn.preprocessing import KBinsDiscretizer, StandardScaler, OneHotEncoder, OrdinalEncoder
from sklearn.compose import ColumnTransformer
from sklearn.decomposition import PCA
from sklearn.model_selection import GridSearchCV, StratifiedKFold, cross_validate
from sklearn.metrics import RocCurveDisplay, auc
from sklearn.ensemble import GradientBoostingClassifier
from sklearn.tree import DecisionTreeClassifier
from sklearn.ensemble import RandomForestClassifier
from sklearn.svm import SVC
from xgboost import XGBClassifier
from sklearn.neighbors import KNeighborsClassifier

# pipeline from imblearn to allow for oversampling!
from imblearn.pipeline import Pipeline
from imblearn.over_sampling import ADASYN, SMOTE

from joblib import dump, load

from hyperopt import hp, tpe, fmin, Trials, space_eval, STATUS_OK
from hyperopt.pyll.base import scope
from hyperopt.early_stop import no_progress_loss

import sys
sys.path.insert(1, '../src/')
from produce_final_dataset import read_all_data
 
########################################################################################
# ======================================================================================
# ======================================================================================
########################################################################################

def calc_intervals_criticality(df_graph_features, n_bins, kind):
    '''
    calculate `n_bins` intervals for the determination of critical residues
    
    kind can be:
    - uniform: use uniform binning;
    - kmeans: use kmeans binning
    
    returns the intervals
    '''
    
    if kind == "uniform":
        
        lb_intervals = sorted(pd.cut(df_graph_features["log_betweenness"], 
                             bins=n_bins, include_lowest=True).unique())

        degree_intervals = sorted(pd.cut(df_graph_features["degree"], 
                                  bins=n_bins, include_lowest=True).unique())
        
    elif kind == "kmeans":
        
        X = df_graph_features[["log_betweenness", "degree"]]

        kbd = KBinsDiscretizer(n_bins=n_bins, encode='ordinal', strategy='kmeans').fit(X)

        lb_intervals = kbd.bin_edges_[0]
        lb_intervals = [pd.Interval(lb_intervals[i], 
                                    lb_intervals[i+1],
                                    closed="both") for i in range(len(lb_intervals)-1)]

        degree_intervals = kbd.bin_edges_[1]
        degree_intervals = [pd.Interval(degree_intervals[i],
                                        degree_intervals[i+1],
                                        closed="both") for i in range(len(degree_intervals)-1)]
        
    return lb_intervals, degree_intervals


########################################################################################
# ======================================================================================
# ======================================================================================
########################################################################################

def classify_critical(row, lb_intervals, degree_intervals):
    
    lb = row["log_betweenness"]
    deg = row["degree"]
    
    lb_lower, *_, lb_upper = lb_intervals
    degree_lower, *_, degree_upper = degree_intervals
    
    if deg in degree_lower and lb in lb_lower:
        ans = "LDLB"
        
    elif deg in degree_upper and lb in lb_upper:
        ans = "HDHB"
        
    elif deg in degree_lower and lb in lb_upper:
        ans = "LDHB"
        
    else: 
        ans = "None"
        
    return ans


########################################################################################
# ======================================================================================
# ======================================================================================
########################################################################################

def calc_criticality(df_graph_features, n_bins, kind):
    '''
    calculate the series with criticality of residues
    
    kind can be:
    - uniform: use uniform binning;
    - kmeans: use kmeans binning
    
    returns the intervals
    '''
    
    lb_intervals, degree_intervals = calc_intervals_criticality(df_graph_features, n_bins, kind)
    
    df_graph_features["criticality"] = df_graph_features.apply(lambda row: 
                                                             classify_critical(row, lb_intervals, degree_intervals),
                                                             axis=1)

    # distribution of criticality
    display(df_graph_features["criticality"].value_counts())
    
    return df_graph_features, lb_intervals, degree_intervals


########################################################################################
# ======================================================================================
# ======================================================================================
########################################################################################

def plot_criticality(df_graph_features, lb_intervals, degree_intervals):
    '''
    plot criticality in degre x log betweenness space
    '''
    
    ax = sns.scatterplot(data = df_graph_features, x="degree", y="log_betweenness", hue="criticality")

    # =============================================

    lb_lower, *_, lb_upper = lb_intervals
    degree_lower, *_, degree_upper = degree_intervals

    plt.axhline(y=lb_lower.left, color="k", ls=":", alpha=0.7)
    plt.axhline(y=lb_lower.right, color="k", ls=":", alpha=0.7)
    plt.axhline(y=lb_upper.left, color="k", ls=":", alpha=0.7)
    plt.axhline(y=lb_upper.right, color="k", ls=":", alpha=0.7)

    plt.axvline(x=degree_lower.left, color="k", ls=":", alpha=0.7)
    plt.axvline(x=degree_lower.right, color="k", ls=":", alpha=0.7)
    plt.axvline(x=degree_upper.left, color="k", ls=":", alpha=0.7)
    plt.axvline(x=degree_upper.right, color="k", ls=":", alpha=0.7)

    # =============================================

    ax.legend(handles=ax.legend_.legendHandles, 
              labels=[t.get_text() for t in ax.legend_.texts],
              title=ax.legend_.get_title().get_text(),
              bbox_to_anchor=(1, 0.5), loc="center left")

    plt.show()

########################################################################################
# ======================================================================================
# ======================================================================================
########################################################################################

def is_pareto_efficient(costs):
    """
    Find the pareto-efficient points
    :param costs: An (n_points, n_costs) array
    :return: A (n_points, ) boolean array, indicating whether each point is Pareto efficient
    """
    is_efficient = np.ones(costs.shape[0], dtype = bool)
    
    for i, c in enumerate(costs):
        
        if is_efficient[i]:
            
            # Keep any point with a greater cost
            is_efficient[is_efficient] = np.any(costs[is_efficient] > c, axis=1)  
            is_efficient[i] = True 
            
    return is_efficient


########################################################################################
# ======================================================================================
# ======================================================================================
########################################################################################

def calc_supercrit(df_graph_features):
    '''
    calculate supercriticality based on the pareto front
    of residues using degree, log-betweenness and closeness
    '''

    supercritical_features = df_graph_features["degree log_betweenness closeness".split()].to_numpy()

    supercritical_mask = is_pareto_efficient(supercritical_features)

    ################################################

    df_graph_features["supercritical"] = False

    df_graph_features.loc[df_graph_features.iloc[np.where(supercritical_mask)].index, "supercritical"] = True
    
    return df_graph_features

########################################################################################
# ======================================================================================
# ======================================================================================
########################################################################################

def plot_pareto_front(df_graph_features, 
                      x_plot="degree", y_plot="log_betweenness", z_plot="closeness", 
                      pairplot=True, plot_3d=True):
    
    if pairplot:
        
        sns.pairplot(df_graph_features["degree log_betweenness closeness supercritical".split()], 
                     hue="supercritical"); 

    if plot_3d:
        
        fig = plt.figure(figsize=(8, 8))
        ax = fig.add_subplot(projection='3d')

        x = df_graph_features[x_plot]
        y = df_graph_features[y_plot]
        z = df_graph_features[z_plot]

        ax.scatter(x, y, z, marker="o", alpha=0.5)

        x_sc = df_graph_features.query("supercritical")[x_plot]
        y_sc = df_graph_features.query("supercritical")[y_plot]
        z_sc = df_graph_features.query("supercritical")[z_plot]

        ax.scatter(x_sc, y_sc, z_sc, marker="x", color="red", s=50, alpha=1)

        ax.plot_trisurf(x_sc, y_sc, z_sc, linewidth=0, color="red")

        ax.set_xlabel(x_plot)
        ax.set_ylabel(y_plot)
        ax.set_zlabel(z_plot)

        plt.show() 
  
    
########################################################################################
# ======================================================================================
# ======================================================================================
########################################################################################

def calc_criticality_super(df_graph_features):
    '''
    collumn indicating the (super)criticality of residues
    '''
    
    # this collumn will also include the supercritical residues
    df_graph_features["criticality_super"] = df_graph_features["criticality"].copy()

    # here, supercritical residues are marked as so, only if their "criticality" was "None"
    # notice how there are some HDHB which are suprecritical!
    supercritical_idxs = df_graph_features.query("supercritical").index

    aux = (df_graph_features.loc[supercritical_idxs, "criticality_super"]
                            .apply(lambda x: "SC" if x == "None" else x))

    df_graph_features.loc[supercritical_idxs, "criticality_super"] = aux
    
    return df_graph_features


# new 11/04/22: overwrite the supercritical ones
def calc_criticality_super_overwrite_sc(df_graph_features):
    '''
    collumn indicating the (super)criticality of residues
    '''
    # this collumn will also include the supercritical residues
    df_graph_features["criticality_super"] = df_graph_features["criticality"].copy()

    df_graph_features.loc[df_graph_features.query("supercritical").index, "criticality_super"] = "SC"
    
    return df_graph_features

########################################################################################
# ======================================================================================
# ======================================================================================
########################################################################################
# ML model functions below
########################################################################################
# ======================================================================================
# ======================================================================================
########################################################################################

def preprocessing_pipeline(X, y, encoding="onehot", pca=False, pct_var_exp=1):
    
    # ==========================================
    # numerical features transformer
    
    features_nums = X.select_dtypes(include=np.number).columns.tolist()

    # always use standard scaler!
    transf_list_nums = [("std_scaler", StandardScaler())]
    
    if pca:
        transf_list_nums.append(('pca', PCA(n_components=pct_var_exp)))
        
    transf_feat_nums = Pipeline(transf_list_nums)

    # ==========================================
    # categorical features transformer
    
    features_cats = X.select_dtypes(exclude=np.number).columns.tolist()
    
    if len(features_cats) != 0:

        if encoding == "onehot":
            cat_encoding = ("onehot", OneHotEncoder())
        elif encoding == "ordinal":
            cat_encoding = ("ordinal", OrdinalEncoder())

        transf_feat_cats = Pipeline([cat_encoding])


        # ==========================================
        # final pipeline

        pre_processor = ColumnTransformer([("transf_num", transf_feat_nums, features_nums), 
                                           ("transf_cat", transf_feat_cats, features_cats)])
      
    # case of only numerical features. columntransformer is not necessary, but let's keep it for generality
    else:
        
        pre_processor = ColumnTransformer([("transf_num", transf_feat_nums, features_nums)])

    return pre_processor
    

def build_pipeline(estmtr, estimators, X, y, 
                   encoding="onehot", pca=False, pct_var_exp=0.9,
                   oversample=False):
    '''
    - estimators: dict of estimators, or single estimator with set hps
    - encoding: "onehot" or "ordinal" or False
    - pca: if True, performs pca, keeping a number of PCS 
    such that eplains "pct_var_exp" percentage of variance
    - oversample: "adasyn" or False
    '''
    
    estimator = estimators[estmtr] if isinstance(estimators, dict) else estimators

    pre_processor = preprocessing_pipeline(X, y, encoding=encoding, pca=pca, pct_var_exp=pct_var_exp)

    # final pipeline
    # preprocess and estimator mandatory. pca and oversample optional.

    pipeline_list = [('pre_process', pre_processor)]

    if oversample:
        pipeline_list.append(('adasyn', ADASYN(sampling_strategy='auto', random_state=42)))

    pipeline_list.append((estmtr, estimator))   

    pipe = Pipeline(pipeline_list)
    
    return pipe
    
########################################################################################
# ======================================================================================
# ======================================================================================
########################################################################################

def modeling_pipeline(estmtr, estimators, param_grids, X, y, 
                      encoding="onehot", pca=False, pct_var_exp=0.9,
                      oversample=False):
    '''
    - estimators: dict of estimators, or single estimator with set hps
    - encoding: "onehot" or "ordinal" or False
    - pca: if True, performs pca, keeping a number of PCS 
    such that eplains "pct_var_exp" percentage of variance
    - oversample: "adasyn" or False
    '''
    
    if encoding not in ["onehot", "ordinal", False]:
        raise ValueError("Encoding strategy not available!")
            
    if oversample not in ["adasyn", False]:
        assert ValueError("Oversampling strategy not available!")
        
    # ==========================================
    
    pca_label = "withPCA" if pca else "noPCA"
                   
    ovsmpl_label = oversample if oversample else "noOverSampling"
    
    encoding_label = encoding if encoding else "onlyNumFeats"
    
    results_file_name = f"../data/grid_cv_results/grid_cv_results_{estmtr}_{encoding_label}_{pca_label}_{ovsmpl_label}.parquet"
    
    # ==========================================
    
    # if grid was already completed before, just read file with results
    try:
        
        results_df_stmtr = pd.read_parquet(results_file_name)
        
        # ==========================================

        results_best = results_df_stmtr.iloc[[0]]
        
        print(f"\nGrid results read for {estmtr}!")

    # run the grid
    except:
    
        print(f"\nStarting grid search for {estmtr}")

        pipe = build_pipeline(estmtr, estimators, X, y, 
                              encoding, pca, pct_var_exp,
                              oversample)
        
        param_grid = param_grids[estmtr]

        # ==========================================
        # GRID SEARCH

        cv = StratifiedKFold(n_splits=10)
        
        grid = GridSearchCV(pipe, param_grid, scoring="roc_auc", cv=cv, verbose=10, n_jobs=-1)

        grid.fit(X, y)

        # ==========================================

        results_df_stmtr = pd.DataFrame(grid.cv_results_).sort_values("rank_test_score")

        # dropping estimator_specific parameters columns
        cols = [x for x in results_df_stmtr.columns if "__" not in x]
        results_df_stmtr = results_df_stmtr[cols].copy()

        # model id columns
        results_df_stmtr["estimator"] = estmtr
        results_df_stmtr["encoding"] = encoding
        results_df_stmtr["pca"] = pca
        results_df_stmtr["oversample"] = oversample
        
        model_file_name = f"../models/intermediate_best_models/best_{estmtr}_{encoding_label}_{pca_label}_{ovsmpl_label}.joblib"
        results_df_stmtr["best_model_path"] = model_file_name

        # ==========================================

        results_best = results_df_stmtr.iloc[[0]]

        # ==========================================

        results_df_stmtr.to_parquet(results_file_name)

        print("\nGrid results exported!")
        
        # ==========================================
        # saving the best model
        
        dump(grid.best_estimator_, model_file_name)

    return results_df_stmtr, results_best

########################################################################################
# ======================================================================================
# ======================================================================================
########################################################################################

def bayes_opt(X, y, estimator, estimator_label, hps_space,
              encoding="onehot",
              pca=False, pct_var_exp=0.9,
              oversample=False,
              n_evals=30):
    '''
    - estimator: string identifying estimator to be used. can be:
        - "gb" for GradientBoostingClassifier;
        - "svm" for SVC;
    - etimator_label: label used to identify the estimator as a full string "Gradient Boosting", for example.
    - n_evals: number of hyperopt samples
    - encoding: "onehot" or "ordinal" or False
    '''
    
    allowed_estimators = "gb dt rf xgb svm knn".split()
    
    if estimator not in allowed_estimators:
        raise ValueError(f"Unvailable estimator!\nMust be one of these: {allowed_estimators}")
        
         
    # ============================================================
    # ============================================================
    
    pre_processor = preprocessing_pipeline(X, y, encoding=encoding, pca=pca, pct_var_exp=pct_var_exp)
    
    pipeline_list = [('pre_process', pre_processor)]

    if oversample:
        pipeline_list.append(('adasyn', ADASYN(sampling_strategy='auto', random_state=42)))
    
    # ============================================================
    # ============================================================
    
    
    def objective(hps_space):
        
        # so that the original list is not appended everytime the function is called
        pipeline_list_func = pipeline_list.copy()
        
        if estimator=="gb": estimator_obj = GradientBoostingClassifier(**hps_space, random_state=42)
            
        elif estimator=="dt": estimator_obj = DecisionTreeClassifier(**hps_space, random_state=42)
            
        elif estimator=="rf": estimator_obj = RandomForestClassifier(**hps_space, random_state=42, n_jobs=-1)
            
        elif estimator=="xgb": estimator_obj = XGBClassifier(**hps_space, random_state=42, 
                                                             use_label_encoder=False, n_jobs=-1,
                                                             eval_metric='logloss')
            
        elif estimator=="svm": estimator_obj = SVC(**hps_space, random_state=42)
            
        elif estimator=="knn": estimator_obj = KNeighborsClassifier(**hps_space, n_jobs=-1)
            
        pipeline_list_func.append((estimator, estimator_obj))

        pipe = Pipeline(pipeline_list_func)
    
        # ==================================
        
        print(f"Testing now: {pipe[-1]}\n")

        splitter = StratifiedKFold(n_splits=10, shuffle=True, random_state=42)

        results = cross_validate(estimator=pipe,
                                 X=X, y=y,
                                 cv=splitter,
                                 scoring="roc_auc",
                                 n_jobs=-1)

        mean_score_cv = results['test_score'].mean()
        var_score_cv = results['test_score'].var()

#         return -mean_f1_cv

        return {'status': STATUS_OK,
                'loss': -mean_score_cv,
                'loss_variance': var_score_cv}

    # ============================================================
    # ============================================================

    trials = Trials()

    best_hps = fmin(objective,
                    space=hps_space,
                    algo=tpe.suggest,
                    max_evals=n_evals,
                    trials=trials,
                    early_stop_fn=no_progress_loss(10),
                    # rstate=np.random.default_rng(42),
                   )

    # ============================================================
    # ============================================================

    best_hps = space_eval(hps_space, best_hps)

    print(f"\nBest hyperparameters:\n{best_hps}")

    # ============================================================
    # ============================================================
    
    if estimator=="gb": model_best = GradientBoostingClassifier(**best_hps, random_state=42)
        
    elif estimator=="dt": model_best = DecisionTreeClassifier(**best_hps, random_state=42)
        
    elif estimator=="rf": model_best = RandomForestClassifier(**best_hps, random_state=42, n_jobs=-1)
        
    elif estimator=="xgb": model_best = XGBClassifier(**best_hps, random_state=42, 
                                                         use_label_encoder=False, n_jobs=-1,
                                                         eval_metric='logloss')
        
    elif estimator=="svm": model_best = SVC(**best_hps, random_state=42)
        
    elif estimator=="knn": model_best = KNeighborsClassifier(**best_hps, n_jobs=-1)
        
    # same pipeline, now with best hps for the last step (model)
    pipeline_list.append((estimator, model_best))
    
    pipe_best = Pipeline(pipeline_list)

    # ============================================================
    # ============================================================
    
    splitter = StratifiedKFold(n_splits=10, shuffle=True, random_state=42)

    results = cross_validate(estimator=pipe_best,
                             X=X, y=y,
                             cv=splitter,
                             scoring="roc_auc",
                             return_train_score=True,
                             n_jobs=-1)

    results_df = pd.DataFrame(results)

    plot_cv_roc_from_estimator(pipe_best, estimator_label, X, y) 

    print(results_df["train_score"].describe())
    print()
    print(results_df["test_score"].describe())
    
########################################################################################
# ======================================================================================
# ======================================================================================
########################################################################################

def string_to_dict(d):
    
    return ast.literal_eval(d)


def get_best_hps(estmtr, results_best):
    
#     return string_to_dict((results_best.query(f"estimator == '{estmtr}'")["params"]
#                                        .squeeze()
#                                        .replace(f"{estmtr}__", "")))

    # new 16/10: using .iloc[0], because we may have more than one "estmtr".
    # since results_best is ordered in descending order by "mean_test_score", 
    # we get the first appearing result, hence the .iloc[0]
    return string_to_dict((results_best.query(f"estimator == '{estmtr}'")["params"]
                                       .iloc[0]
                                       .replace(f"{estmtr}__", "")))

def build_pipes_best_hps(estmtr, estimators_best_hps, X, y, results_best):

    encoding, pca, oversample = results_best.query(f"estimator == '{estmtr}'")[["encoding",
                                                                                "pca",
                                                                                "oversample"]].iloc[0]

    pipe = build_pipeline(estmtr, estimators_best_hps[estmtr], X, y, 
                          encoding=encoding, pca=pca, oversample=oversample)
    
    return pipe

########################################################################################
# ======================================================================================
# ======================================================================================
########################################################################################

def get_best_pipelines(X, y, results_best):

    # new 16/10/22: since we're saving the best models (and their path is in results_best), 
    # it's much simpler (and safe) to do as below!

    # this gets the joblib file for each one of the best models, after the HPs optimization
    # since results_best is ordered in descending order by "mean_test_score", 
    # we get the first appearing result, hence the idx[0]!
    best_models = {estmtr : results_best.loc[idx[0], "best_model_path"] 
                   for estmtr, idx in results_best.groupby("estimator").groups.items()}

    # read the model objects
    pipelines_best_hps = {estmtr : load(path) for estmtr, path in best_models.items()}

    return pipelines_best_hps

########################################################################################
# ======================================================================================
# ======================================================================================
########################################################################################

def plot_cv_roc_from_estimator(estimator, estimator_label, X, y):
    '''
    pass only estimator with best hyper params, not trained models
    '''

    cv = StratifiedKFold(n_splits=10)

    tprs = []
    aucs = []
    mean_fpr = np.linspace(0, 1, 100)

    fig, ax = plt.subplots(figsize=(8, 6))

    for i, (train, test) in enumerate(cv.split(X, y)):
        
        model = estimator.fit(X.iloc[train], y.iloc[train])
        
        viz = RocCurveDisplay.from_estimator(model, 
                                             X.iloc[test], y.iloc[test],
                                             alpha=0.3, lw=1, 
                                             label=False,
                                             ax=ax)

        interp_tpr = np.interp(mean_fpr, viz.fpr, viz.tpr)
        interp_tpr[0] = 0.0

        tprs.append(interp_tpr)
        aucs.append(viz.roc_auc)

    # diagonal line
    ax.plot([0, 1], [0, 1], linestyle="--", lw=2, color="k", label="Random", alpha=0.8)

    mean_tpr = np.mean(tprs, axis=0)
    mean_tpr[-1] = 1.0
    mean_auc = auc(mean_fpr, mean_tpr)
    std_auc = np.std(aucs)
    ax.plot(mean_fpr, mean_tpr, color="b",
            label=r"Mean ROC (AUC = %0.2f $\pm$ %0.2f)" % (mean_auc, std_auc),
            lw=2, alpha=0.8)

    std_tpr = np.std(tprs, axis=0)
    tprs_upper = np.minimum(mean_tpr + std_tpr, 1)
    tprs_lower = np.maximum(mean_tpr - std_tpr, 0)
    ax.fill_between(mean_fpr, tprs_lower, tprs_upper, color="grey",
                    alpha=0.2, label=r"$\pm$ 1 std. dev.")

    ax.set(xlim=[-0.05, 1.05],
           ylim=[-0.05, 1.05])

    plt.title(f"ROC curve - {estimator_label}", fontsize=18)

    plt.gcf()
    handles, labels = plt.gca().get_legend_handles_labels()
    handles = handles[10:]
    labels = labels[10:]
    by_label = dict(zip(labels, handles))
    plt.legend(by_label.values(), by_label.keys(), loc="lower right")
    plt.show()
    
########################################################################################
# ======================================================================================
# ======================================================================================
########################################################################################

def create_roc_curve(X, y, estimator, ax, label, color):
    cv = StratifiedKFold(n_splits=10)

    tprs = []
    aucs = []
    mean_fpr = np.linspace(0, 1, 100)

    for i, (train, test) in enumerate(cv.split(X, y)):

        model = estimator.fit(X.iloc[train], y.iloc[train])

        viz = RocCurveDisplay.from_estimator(model, 
                                             X.iloc[test], y.iloc[test],
                                             alpha=0.3, lw=1, 
                                             label=False)
        plt.close()
    
        interp_tpr = np.interp(mean_fpr, viz.fpr, viz.tpr)
        interp_tpr[0] = 0.0

        tprs.append(interp_tpr)
        aucs.append(viz.roc_auc)
    

    mean_tpr = np.mean(tprs, axis=0)
    mean_tpr[-1] = 1.0
    mean_auc = auc(mean_fpr, mean_tpr)
    std_auc = np.std(aucs)
    ax.plot(mean_fpr, mean_tpr,
            label="%ls: %0.2f $\pm$ %0.2f" % (label.upper(), mean_auc, std_auc),
            color=color,
            lw=2, alpha=0.8)

    std_tpr = np.std(tprs, axis=0)

    ax.set(xlim=[-0.05, 1.05],
           ylim=[-0.05, 1.05])
    
########################################################################################
# ======================================================================================
# ======================================================================================
########################################################################################

# bringing residue name (3 letter code)
def get_res_name(row):
    
    res1, res2 = row[["Res1", "Res2"]]
    
    if res1 == res2:
        ans = res1
    else:
        if str(res1) == "nan":
            ans = res2
        elif str(res2) == "nan":
            ans = res1
        else:
            raise ValueError("error, at least one must be not null!")
    
    return ans

def get_name(row):
    
    node = row["node"]
    residue = row["residue"]
    
    # format of return: capitalized residue name (3-lettercode), then the residue number
    # for instance, Ala1 (Alanine, first node)
    return f"{residue.capitalize()}{node}"

def graph_features_with_res_name_code(df_graph_features, prot, str_uni, str_inter):
    
    # this table has the residues name (3 letter code)
    graph_data_path = os.path.join("../data/igraph/", 'out', f"{prot}_graph_data_{str_uni}{str_inter}.csv")
    
    df_graph_data = pd.read_csv(graph_data_path)

    # merge to bring the residue name
    df_graph_features = (df_graph_features.merge(df_graph_data[["Pos1", "Res1"]], 
                                                 left_on="node", right_on="Pos1", how="left")
                                          .drop_duplicates()
                                          .merge(df_graph_data[["Pos2", "Res2"]], 
                                                 left_on="node", right_on="Pos2", how="left")
                                          .drop_duplicates())

    # function to make sure the correct residue name is captured
    df_graph_features["residue"] = df_graph_features.apply(get_res_name, axis=1)

    return df_graph_features

def get_nodes_with_zero_betweenness_formated(df_graph_features, prot, str_uni, str_inter):

    df_graph_features = graph_features_with_res_name_code(df_graph_features, 
                                                          prot, 
                                                          str_uni,
                                                          str_inter)

    # list of nodes with zero betweenness, elements in the desired format
    nodes_zero_betw = (df_graph_features.query("betweenness == 0")[["node", "residue"]]
                                        .sort_values("node")
                                        .apply(get_name, axis=1)
                                        .values.tolist())
    
    print(f"There are {len(nodes_zero_betw)} nodes with zero betweenness in the RIN!\n")

    return nodes_zero_betw

def get_disconnected_nodes_formated(df_graph_features, disconnected_res, prot, str_uni, str_inter):

    df_graph_features = graph_features_with_res_name_code(df_graph_features, 
                                                          prot, 
                                                          str_uni,
                                                          str_inter)
    
    # list of disconnected nodes, elements in the desired format
    nodes_disconnected = (df_graph_features.query(f"node in {disconnected_res}")
                                           .sort_values("node")
                                           .apply(get_name, axis=1)
                                           .values.tolist())
    
    print(f"There are {len(nodes_disconnected)} disconnected nodes in the RIN!\n")

    return nodes_disconnected

########################################################################################
# ======================================================================================
# ======================================================================================
########################################################################################

def plot_criticality_pretty(df, results_figs_path):
    
    g = sns.JointGrid(data=df, x='degree', y='log_betweenness', hue='mutation_observed', height=10,ratio=5)
    color_dict = {0: to_rgba('black', 0.1),
                  1: to_rgba('red', 1)}

    kwargs  =   {'edgecolor':"k", # for edge color
                 'linewidth':0.45, # line width of spot
                 'linestyle':'-', # line style of spot
                }
    g.plot_joint(sns.scatterplot, palette=color_dict, lw=10, **kwargs)
    g.plot_marginals(sns.histplot, kde=False, palette=color_dict, multiple='stack')

    scatter_ax = g.fig.get_axes()[0]
    scatter_ax.set_xlabel("Degree", fontsize=16)
    scatter_ax.set_ylabel("Log Betwenness", fontsize=16)


    def set_mutation_legend(ax):
        """
        Calls legend, and sets all the legend colors opacity to 100%.
        Returns the legend handle.
        """
        leg = ax.legend(bbox_to_anchor=(1,0.1),loc='lower right', fontsize=12, title='Mutation Database Presence', ncol=2)
        for lh in leg.legendHandles:
            fc_arr = lh.get_fc().copy()
            fc_arr[:, -1] = 1
            lh.set_fc(fc_arr)

        return leg

    leg = set_mutation_legend(scatter_ax)


    def classify_critical(row, lb_intervals, degree_intervals):

        lb = row["log_betweenness"]
        deg = row["degree"]

        lb_lower, *_, lb_upper = lb_intervals
        degree_lower, *_, degree_upper = degree_intervals

        if deg in degree_lower and lb in lb_lower:
            ans = "LDLB"

        elif deg in degree_upper and lb in lb_upper:
            ans = "HDHB"

        elif deg in degree_lower and lb in lb_upper:
            ans = "LDHB"

        else: 
            ans = "None"

        return ans

    n_bins = 4

    lb_intervals = sorted(pd.cut(df["log_betweenness"], 
                                 bins=n_bins, include_lowest=True).unique())

    degree_intervals = sorted(pd.cut(df["degree"], 
                                     bins=n_bins, include_lowest=True).unique())
    lb_lower, *_, lb_upper = lb_intervals
    degree_lower, *_, degree_upper = degree_intervals
    x1, y1 = degree_lower.left-0.15, lb_lower.left
    x2, y2 = degree_lower.left-0.15, lb_upper.left
    x3, y3 = degree_upper.left, lb_upper.left
    r1 = scatter_ax.add_patch(Rectangle((x1, y1), degree_lower.right - x1, lb_lower.right-y1, facecolor="lightskyblue", alpha=0.3, edgecolor="black"))
    r2 = scatter_ax.add_patch(Rectangle((x2, y2), degree_lower.right - x2, lb_upper.right-y2, facecolor="khaki", alpha=0.3, edgecolor="black"))
    r3 = scatter_ax.add_patch(Rectangle((x3, y3), degree_upper.right - x3 + 0.1, lb_upper.right-y3, facecolor="rosybrown", alpha=0.3, edgecolor="black"))


    criticity_ax = scatter_ax.twinx()
    criticity_legend = plt.legend([r1, r2, r3], ["LDLB", "LDHB", "HDHB"], loc='lower right', fontsize=12, title='Node Criticity Analysis', ncol=3, frameon=True)
    for lh in criticity_legend.legendHandles: 
        lh.set_alpha(1)
    criticity_ax.add_artist(criticity_legend)
    criticity_ax.get_yaxis().set_visible(False)

    plt.tight_layout()
    plt.savefig(f"{results_figs_path}/criticality_figure.png", format="png")
    plt.show()
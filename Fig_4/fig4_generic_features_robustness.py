import pandas as pd
import numpy as np
from scipy.stats import spearmanr
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib as mpl


def calc_percent_significant_asvs(model,data,data_):
    results = pd.read_csv("rev_fig/consistency_of_model_between_datasets/all_corrected_fp_tp.csv",index_col=0)
    num_significant = results[f"TP-{model}"][data]
    total_data = len(pd.read_csv(f"ancom_bc2/data/{data_}/data_for_deseq_{data_}.csv",index_col=0))
    return num_significant/total_data

def calc_rsp(model,data,beta):
    results = pd.read_csv("rev_fig/consistency_of_model_between_datasets/all_corrected_fp_tp.csv", index_col=0)
    RP = results[f"TP-{model}"][data]
    SP = results[f"FP-{model}"][data]
    RSP = (beta * RP - SP) / (beta * RP + SP)
    return RSP


if __name__ == '__main__':
    SIZE = 20
    mpl.rc('font', family='Times New Roman')
    GENERIC = ["log(sample size)","Read depth variation","log(read depth)","Sparsity","% below prev. cut-off","Richness"]
    LIST_MODELS = ["miMic","miMic relative","LINDA","LINDA-C","DeSeq","DeSeq-C","ANCOM","ANCOM-C","ANCOM-BC2","ANCOM-BC2-C","ALDEx Welch","ALDEx Welch-C","ALDEx Wilcoxon","ALDEx Wilcoxon-C","LEfSe"]
    LIST_DATASETS_16S = ["Cirrhosis_Knight_Lab","ERP020401","ERP021216","IBD","he","jacob","MF","ok_30","ok_94","PRJNA353587","PRJNA419097","WB"]
    LIST_DATASETS_WGS = ["Cirrhosis","data1_metab","data2_metab","data3_metab","data5_metab","data6_metab","T2D","obesity"]
    WGS = False
    RSP = False
    if WGS:
        LIST_DATASETS = LIST_DATASETS_WGS
    else:
        LIST_DATASETS = LIST_DATASETS_16S
    generic_features = pd.DataFrame(index=LIST_DATASETS,columns=GENERIC+LIST_MODELS)
    for DATA in LIST_DATASETS:
        df = pd.read_csv(f"ancom_bc2/data/{DATA}/data_for_deseq_{DATA}.csv",index_col=0)

        # Calculate features
        # RICHNESS
        richness = (df.astype(bool).sum(axis=0)).mean()
        log_sample_size = np.log10(len(df))
        read_depth_var = df.sum(axis=1).var()
        log_read_depth = np.log10(df.sum(axis=1).median())

        # sparsity
        total_entries = df.shape[0]*df.shape[1]
        sparsity = ((df == 0.0).sum().sum())/(df.shape[0]*df.shape[1])

        # percent of ASVs under 10% prevalence
        # Step 1: Calculate prevalence of each feature
        prevalence = (df.astype(bool).sum(axis=0) / len(df)) * 100
        # Step 2: Count number of features with prevalence under 10%
        under_10_percent = (prevalence < 10).sum()
        percent_under_10_percent = (under_10_percent / len(df.columns)) * 100

        # Update pandas
        generic_features["log(sample size)"][DATA] = log_sample_size
        generic_features["Read depth variation"][DATA] = read_depth_var
        generic_features["log(read depth)"][DATA] = log_read_depth
        generic_features["Sparsity"][DATA] = sparsity
        generic_features["% below prev. cut-off"][DATA] = percent_under_10_percent
        generic_features["Richness"][DATA] = richness

    # Build model data
    if WGS:
        generic_features = generic_features.rename(index={'Cirrhosis': 'Cirrhosis-POP', "obesity":"Obesity","data1_metab":"WGS-1","data2_metab":"WGS-2","data3_metab":"WGS-3","data5_metab":"WGS-4","data6_metab":"WGS-5"})
    else:
        generic_features = generic_features.rename(
            index={'Cirrhosis_Knight_Lab': 'Cirrhosis','he':'He','jacob':'Jacob',"ok_30":"OK30","ok_94":"OK94"} )
    for model in LIST_MODELS:
        for (data,data_) in zip(list(generic_features.index),LIST_DATASETS):
            if RSP:
                generic_features[model][data] = calc_rsp(model,data,1)
            else:
                generic_features[model][data] = calc_percent_significant_asvs(model,data,data_)

    # Calculate correlation
    generic_features = generic_features.fillna(0.0)
    corrs_df = pd.DataFrame(index=LIST_MODELS, columns=GENERIC)
    p_values_df = pd.DataFrame(index=LIST_MODELS, columns=GENERIC)

    for model in LIST_MODELS:
        for g in GENERIC:
            scc, p = spearmanr(generic_features[g], generic_features[model])
            corrs_df[g][model] = scc
            p_values_df[g][model] = p
    corrs_df = corrs_df.fillna(0.0)
    corrs_df = corrs_df.T
    # Plot heatmap
    plt.figure(figsize=(11, 5))
    heatmap = sns.heatmap(corrs_df.astype(float), cmap='coolwarm',fmt=".2f", vmin=-1.0, vmax=1.0)#annot=True,

    # Add asterisks for significant correlations
    for i in range(len(LIST_MODELS)):
        for j in range(len(GENERIC)):
            if p_values_df.iloc[i, j] < 0.05:
                plt.text(i + 0.5, j + 0.5, "*", ha='center', va='center', color='black', fontsize=12)

    plt.ylabel('Generic features',fontsize=SIZE)
    plt.xlabel('Models',fontsize=SIZE)
    plt.xticks(fontsize=SIZE)
    plt.yticks(fontsize=SIZE)
    if RSP:
        if WGS:
            plt.title(
                "Dataset characteristics (WGS) associated with\nRSP score",fontsize=SIZE,fontweight="bold")
        else:
            plt.title(
                "Dataset characteristics (16S) associated with\nRSP score",fontsize=SIZE,fontweight="bold")
    else:
        if WGS:
            plt.title(
                "Dataset characteristics (WGS) associated with\npercentage of significant amplicon sequence variants",fontsize=SIZE,fontweight="bold")
        else:
            plt.title("Dataset characteristics (16S) associated with\npercentage of significant amplicon sequence variants",fontsize=SIZE,fontweight="bold")


    # Add colors to xticks labels
    ax = plt.gca()
    ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha='right', fontsize=20)

    # Set colors for xticks labels
    LIST_COLORS = [ "hotpink", "mediumorchid","lightskyblue", "cornflowerblue", "salmon", "tomato", "mediumaquamarine",
                   "seagreen","green", "darkgreen", "grey",
                   "dimgrey", "lightslategray", "slategray",  'gold']
    for ticklabel, color in zip(ax.get_xticklabels(), LIST_COLORS):
        ticklabel.set_color(color)

        # Increase font size of colorbar
    cbar = heatmap.collections[0].colorbar
    cbar.ax.tick_params(labelsize=20)

    plt.tight_layout()

    if RSP:
        if WGS:
            plt.savefig("rev_fig/only_significant/heatmap_generic_rsp_wgs.png")
        else:
            plt.savefig("rev_fig/only_significant/heatmap_generic_rsp_16s.png")
    else:
        if WGS:
            plt.savefig("rev_fig/only_significant/heatmap_generic_percent_tp_wgs.png")
        else:
            plt.savefig("rev_fig/only_significant/heatmap_generic_percent_tp_16s.png")

    plt.show()




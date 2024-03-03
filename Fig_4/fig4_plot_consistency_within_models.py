import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl

def plot_consistency_model(merged,model,color,merge,avg=False):
    SIZE = 20
    if not avg:
        significant_per_model = merged[merged[model] == True]
        num_s = len(significant_per_model)
        counts = significant_per_model['Total'].value_counts()
        counts = counts.to_frame('Counts')
        # Create a DataFrame with the desired index values
        if merge:
            index_values = [1, 2, 3, 4, 5, 6, 7]
        else:
            index_values = [1, 2, 3, 4, 5,6,7,8,9,10,11,12,13,14]
        counts_df = pd.DataFrame(index=index_values)
        counts_df = pd.concat([counts_df, counts], axis=1)
        counts_df = counts_df.fillna(0.0)
        counts_df = (counts_df / num_s) * 100
        counts_df.plot(kind="bar", color=color, rot=0, figsize=(5, 3))
    else:
        counts_df = merged
        counts_df.plot(kind="bar", color=color, rot=0, figsize=(5, 5),edgecolor="black")
    plt.title(model, fontsize=SIZE, fontweight="bold")
    plt.xlabel("No. tools that called\nthe feature significant",fontsize=SIZE)
    plt.ylabel("Percent total\nsignificant features",fontsize=SIZE)
    plt.xticks(fontsize=SIZE)
    plt.yticks(fontsize=SIZE)
    plt.ylim([0,80])
    #plt.yscale("log")
    if model != "miMic":
        if not avg:
            plt.legend([f"Total hits = {num_s}"], fontsize=SIZE)
        else:
            plt.legend([f"Total hits = 509"], fontsize=SIZE)
    else:
        if not avg:
            plt.legend([f"Total hits = {num_s+24}"], fontsize=SIZE)
        else:
            plt.legend([f"Total hits = 390"], fontsize=SIZE)
    plt.tight_layout()
    if merge:
        if avg:
            plt.savefig(f"rev_fig/consistent_model/merge/avg_{model}.png")
        else:
            plt.savefig(f"rev_fig/consistent_model/merge/{model}.png")
    else:
        plt.savefig(f"rev_fig/consistent_model/no_merge/{model}.png")
    plt.show()
    if not avg:
        return counts_df

def plot_paired_consistency_model(avg_to_plot,mimic,c1,c2):
    SIZE = 25
    to_plot = pd.DataFrame(index=mimic.index,columns=["miMic\n(Total hits = 390)","Average over all models\nTotal hits = 509"])
    to_plot["miMic\n(Total hits = 390)"] = mimic
    to_plot["Average over all models\nTotal hits = 509"] = avg_to_plot
    to_plot.plot(kind="bar",color=[c2,c1],edgecolor="black",figsize=(10, 5),rot=0)
    plt.xlabel("No. tools that called\nthe feature significant", fontsize=SIZE)
    plt.ylabel("Percent total\nsignificant features", fontsize=SIZE)
    plt.xticks(fontsize=SIZE)
    plt.yticks(fontsize=SIZE)
    plt.legend(fontsize=SIZE)
    plt.ylim([0, 80])
    plt.tight_layout()
    plt.savefig("rev_fig/consistent_model/merge/avg_mimic_comp.png")
    plt.show()
    c=0

def add_missing(df,LIST_MODELS):
    for model in LIST_MODELS:
        if model not in list(df.columns):
            df[model] = False
    return df


def build_combined_df(all_models,LIST_MODELS):
    all_models = add_missing(all_models,LIST_MODELS)
    all_models["b-LINDA"] = all_models["LINDA"] | all_models["LINDA-C"]
    all_models["b-miMic"] = all_models["miMic"] | all_models["miMic relative"]
    all_models["b-ANCOM"] = all_models["ANCOM"] | all_models["ANCOM-C"]
    all_models["b-DeSeq"] = all_models["DeSeq"] | all_models["DeSeq-C"]
    all_models["b-ANCOM-BC2"] = all_models["ANCOM-BC2"] | all_models["ANCOM-BC2-C"]
    all_models["b-ALDEx"] = all_models["ALDEx Welch"] | all_models["ALDEx Welch-C"] | all_models["ALDEx Wilcoxon"] | \
                            all_models["ALDEx Wilcoxon-C"]
    return all_models

def load_p_vals(model,data):
    if model == "LINDA":
        df = pd.read_csv(f"linda_results/LINDA_{data}.csv",index_col=0)
        p = df["output.Tag.pvalue"]
        p.index = [i.replace("..",";").replace(".",";") for i in p.index]
    elif model == "LINDA-C":
        df = pd.read_csv(f"linda_results/LINDA_{data}.csv", index_col=0)
        p = df["output.Tag.padj"]
        p.index = [i.replace("..", ";").replace(".", ";") for i in p.index]
    elif model == "DeSeq":
        p = pd.read_csv(f"ancom_bc2/data/IBD/deseq_p_vals.csv", index_col=0)
        p.index = [i.replace("; ", ";") for i in p.index]
        p.index = [i.replace("_p", ";p").replace("_c", ";c").replace("_o", ";o").replace("_f", ";f").replace("_g",
                                                                                                             ";g").replace(
            "_s", ";s") for i in p.index]
        list_of_lists = [(i.replace("k__", "").replace("p__", "").replace("c__", "").replace("o__", "").replace("f__",
                                                                                                                "").replace(
            "g__", "").replace("s__", "")).split(";") for i in p.index]
        list_of_lists_without_empty_strings = [[item for item in sublist if item != ""] for sublist in list_of_lists]
        list_leaves = list_of_lists_without_empty_strings
        list_lens = [len(i) for i in list_leaves]
        otu_train_cols = list()
        for i, j in zip(list_leaves, list_lens):
            if j == 1:
                updated = "k__" + i[0]
            elif j == 2:
                updated = "k__" + i[0] + ";" + "p__" + i[1]
            elif j == 3:
                updated = "k__" + i[0] + ";" + "p__" + i[1] + ";" + "c__" + \
                          i[2]
            elif j == 4:
                updated = "k__" + i[0] + ";" + "p__" + i[1] + ";" + "c__" + i[
                    2] + ";" + "o__" + i[3]
            elif j == 5:
                updated = "k__" + i[0] + ";" + "p__" + i[1] + ";" + "c__" + i[
                    2] + ";" + "o__" + i[3] + ";" + "f__" + i[4]
            elif j == 6:
                updated = "k__" + i[0] + ";" + "p__" + i[1] + ";" + "c__" + i[
                    2] + ";" + "o__" + i[3] + ";" + "f__" + i[4] + ";" + "g__" + \
                          i[5]
            elif j == 7:
                updated = "k__" + i[0] + ";" + "p__" + i[1] + ";" + "c__" + i[
                    2] + ";" + "o__" + i[3] + ";" + "f__" + i[4] + ";" + "g__" + i[
                              5] + ";" + "s__" + i[6]
            otu_train_cols.append(updated)
        p.index = otu_train_cols
        p = p.groupby(p.index).first()

    elif model == "DeSeq-C":
        p = pd.read_csv(f"ancom_bc2/data/{data}/deseq_p_vals_corrected.csv", index_col=0)
        p.index = [i.replace("; ", ";") for i in p.index]
        p.index = [i.replace("_p", ";p").replace("_c", ";c").replace("_o", ";o").replace("_f", ";f").replace("_g",
                                                                                                             ";g").replace(
            "_s", ";s") for i in p.index]
        list_of_lists = [(i.replace("k__", "").replace("p__", "").replace("c__", "").replace("o__", "").replace("f__",
                                                                                                                "").replace(
            "g__", "").replace("s__", "")).split(";") for i in p.index]
        list_of_lists_without_empty_strings = [[item for item in sublist if item != ""] for sublist in list_of_lists]
        list_leaves = list_of_lists_without_empty_strings
        list_lens = [len(i) for i in list_leaves]
        otu_train_cols = list()
        for i, j in zip(list_leaves, list_lens):
            if j == 1:
                updated = "k__" + i[0]
            elif j == 2:
                updated = "k__" + i[0] + ";" + "p__" + i[1]
            elif j == 3:
                updated = "k__" + i[0] + ";" + "p__" + i[1] + ";" + "c__" + \
                          i[2]
            elif j == 4:
                updated = "k__" + i[0] + ";" + "p__" + i[1] + ";" + "c__" + i[
                    2] + ";" + "o__" + i[3]
            elif j == 5:
                updated = "k__" + i[0] + ";" + "p__" + i[1] + ";" + "c__" + i[
                    2] + ";" + "o__" + i[3] + ";" + "f__" + i[4]
            elif j == 6:
                updated = "k__" + i[0] + ";" + "p__" + i[1] + ";" + "c__" + i[
                    2] + ";" + "o__" + i[3] + ";" + "f__" + i[4] + ";" + "g__" + \
                          i[5]
            elif j == 7:
                updated = "k__" + i[0] + ";" + "p__" + i[1] + ";" + "c__" + i[
                    2] + ";" + "o__" + i[3] + ";" + "f__" + i[4] + ";" + "g__" + i[
                              5] + ";" + "s__" + i[6]
            otu_train_cols.append(updated)
        p.index = otu_train_cols
        p = p.groupby(p.index).first()
    elif model == "ANCOM":
        p = pd.read_csv(f"ancom_bc2/data/{data}/ancom_p_vals.csv", index_col=0)
        p = p[p["Reject null hypothesis"] == True]
        p["p"] = p["Reject null hypothesis"]
        p = p["p"]
        p.index = [i.replace("; ", ";") for i in p.index]
    elif model == "ANCOM-C":
        p = pd.read_csv(f"ancom_bc2/data/{data}/ancom_p_vals_corrected.csv", index_col=0)
        p = p[p["Reject null hypothesis"] == True]
        p["p"] = p["Reject null hypothesis"]
        p = p["p"]
        p.index = [i.replace("; ", ";") for i in p.index]
    elif model == "ALDEx Welch":
        p = pd.read_csv(f"ancom_bc2/data/{data}/aldex_real_holm.csv", index_col=0)[f"we.ep"]
        p.index = [i.replace("..", ";").replace(".", ";") for i in p.index]
        p.index = [i.replace("_p", ";p").replace("_c", ";c").replace("_o", ";o").replace("_f", ";f").replace("_g",
                                                                                                             ";g").replace(
            "_s", ";s") for i in p.index]
        if data == "MF":
            p.index = [i.split("t__")[0] for i in p.index]
        list_of_lists = [(i.replace("k__", "").replace("p__", "").replace("c__", "").replace("o__", "").replace("f__",
                                                                                                                "").replace(
            "g__", "").replace("s__", "")).split(";") for i in p.index]
        list_of_lists_without_empty_strings = [[item for item in sublist if item != ""] for sublist in list_of_lists]
        list_leaves = list_of_lists_without_empty_strings
        list_lens = [len(i) for i in list_leaves]
        otu_train_cols = list()
        for i, j in zip(list_leaves, list_lens):
            if j == 1:
                updated = "k__" + i[0]
            elif j == 2:
                updated = "k__" + i[0] + ";" + "p__" + i[1]
            elif j == 3:
                updated = "k__" + i[0] + ";" + "p__" + i[1] + ";" + "c__" + \
                          i[2]
            elif j == 4:
                updated = "k__" + i[0] + ";" + "p__" + i[1] + ";" + "c__" + i[
                    2] + ";" + "o__" + i[3]
            elif j == 5:
                updated = "k__" + i[0] + ";" + "p__" + i[1] + ";" + "c__" + i[
                    2] + ";" + "o__" + i[3] + ";" + "f__" + i[4]
            elif j == 6:
                updated = "k__" + i[0] + ";" + "p__" + i[1] + ";" + "c__" + i[
                    2] + ";" + "o__" + i[3] + ";" + "f__" + i[4] + ";" + "g__" + \
                          i[5]
            elif j == 7:
                updated = "k__" + i[0] + ";" + "p__" + i[1] + ";" + "c__" + i[
                    2] + ";" + "o__" + i[3] + ";" + "f__" + i[4] + ";" + "g__" + i[
                              5] + ";" + "s__" + i[6]
            elif j == 8:
                updated = "k__" + i[0] + ";" + "p__" + i[1] + ";" + "c__" + i[
                    2] + ";" + "o__" + i[3] + ";" + "f__" + i[4] + ";" + "g__" + i[
                              5] + ";" + "s__" + i[6]
            elif j == 9:
                updated = "k__" + i[0] + ";" + "p__" + i[1] + ";" + "c__" + i[
                    2] + ";" + "o__" + i[3] + ";" + "f__" + i[4] + ";" + "g__" + i[
                              5] + ";" + "s__" + i[6]
            elif j == 10:
                updated = "k__" + i[0] + ";" + "p__" + i[1] + ";" + "c__" + i[
                    2] + ";" + "o__" + i[3] + ";" + "f__" + i[4] + ";" + "g__" + i[
                              5] + ";" + "s__" + i[6]
            elif j == 11:
                updated = "k__" + i[0] + ";" + "p__" + i[1] + ";" + "c__" + i[
                    2] + ";" + "o__" + i[3] + ";" + "f__" + i[4] + ";" + "g__" + i[
                              5] + ";" + "s__" + i[6]
            elif j == 12:
                updated = "k__" + i[0] + ";" + "p__" + i[1] + ";" + "c__" + i[
                    2] + ";" + "o__" + i[3] + ";" + "f__" + i[4] + ";" + "g__" + i[
                              5] + ";" + "s__" + i[6]
            otu_train_cols.append(updated)
        p.index = otu_train_cols
        p = p.groupby(p.index).first()
    elif model == "ALDEx Welch-C":
        p = pd.read_csv(f"ancom_bc2/data/{data}/aldex_real_holm.csv",index_col=0)[f"we.eBH"]
        p.index = [i.replace("..", ";").replace(".", ";") for i in p.index]
        p.index = [i.replace("_p", ";p").replace("_c", ";c").replace("_o", ";o").replace("_f", ";f").replace("_g",
                                                                                                             ";g").replace(
            "_s", ";s") for i in p.index]
        if data == "MF":
            p.index = [i.split("t__")[0] for i in p.index]
        list_of_lists = [(i.replace("k__", "").replace("p__", "").replace("c__", "").replace("o__", "").replace("f__",
                                                                                                                "").replace(
            "g__", "").replace("s__", "")).split(";") for i in p.index]
        list_of_lists_without_empty_strings = [[item for item in sublist if item != ""] for sublist in list_of_lists]
        list_leaves = list_of_lists_without_empty_strings
        list_lens = [len(i) for i in list_leaves]
        otu_train_cols = list()
        for i, j in zip(list_leaves, list_lens):
            if j == 1:
                updated = "k__" + i[0]
            elif j == 2:
                updated = "k__" + i[0] + ";" + "p__" + i[1]
            elif j == 3:
                updated = "k__" + i[0] + ";" + "p__" + i[1] + ";" + "c__" + \
                          i[2]
            elif j == 4:
                updated = "k__" + i[0] + ";" + "p__" + i[1] + ";" + "c__" + i[
                    2] + ";" + "o__" + i[3]
            elif j == 5:
                updated = "k__" + i[0] + ";" + "p__" + i[1] + ";" + "c__" + i[
                    2] + ";" + "o__" + i[3] + ";" + "f__" + i[4]
            elif j == 6:
                updated = "k__" + i[0] + ";" + "p__" + i[1] + ";" + "c__" + i[
                    2] + ";" + "o__" + i[3] + ";" + "f__" + i[4] + ";" + "g__" + \
                          i[5]
            elif j == 7:
                updated = "k__" + i[0] + ";" + "p__" + i[1] + ";" + "c__" + i[
                    2] + ";" + "o__" + i[3] + ";" + "f__" + i[4] + ";" + "g__" + i[
                              5] + ";" + "s__" + i[6]
            elif j == 8:
                updated = "k__" + i[0] + ";" + "p__" + i[1] + ";" + "c__" + i[
                    2] + ";" + "o__" + i[3] + ";" + "f__" + i[4] + ";" + "g__" + i[
                              5] + ";" + "s__" + i[6]
            elif j == 9:
                updated = "k__" + i[0] + ";" + "p__" + i[1] + ";" + "c__" + i[
                    2] + ";" + "o__" + i[3] + ";" + "f__" + i[4] + ";" + "g__" + i[
                              5] + ";" + "s__" + i[6]
            elif j == 10:
                updated = "k__" + i[0] + ";" + "p__" + i[1] + ";" + "c__" + i[
                    2] + ";" + "o__" + i[3] + ";" + "f__" + i[4] + ";" + "g__" + i[
                              5] + ";" + "s__" + i[6]
            elif j == 11:
                updated = "k__" + i[0] + ";" + "p__" + i[1] + ";" + "c__" + i[
                    2] + ";" + "o__" + i[3] + ";" + "f__" + i[4] + ";" + "g__" + i[
                              5] + ";" + "s__" + i[6]
            elif j == 12:
                updated = "k__" + i[0] + ";" + "p__" + i[1] + ";" + "c__" + i[
                    2] + ";" + "o__" + i[3] + ";" + "f__" + i[4] + ";" + "g__" + i[
                              5] + ";" + "s__" + i[6]
            otu_train_cols.append(updated)
        p.index = otu_train_cols
        p = p.groupby(p.index).first()
    elif model == "ALDEx Wilcoxon":
        p = pd.read_csv(f"ancom_bc2/data/{data}/aldex_real_holm.csv", index_col=0)[f"wi.ep"]
        p.index = [i.replace("..", ";").replace(".", ";") for i in p.index]
        p.index = [i.replace("_p", ";p").replace("_c", ";c").replace("_o", ";o").replace("_f", ";f").replace("_g",
                                                                                                             ";g").replace(
            "_s", ";s") for i in p.index]
        if data == "MF":
            p.index = [i.split("t__")[0] for i in p.index]
        list_of_lists = [(i.replace("k__", "").replace("p__", "").replace("c__", "").replace("o__", "").replace("f__",
                                                                                                                "").replace(
            "g__", "").replace("s__", "")).split(";") for i in p.index]
        list_of_lists_without_empty_strings = [[item for item in sublist if item != ""] for sublist in list_of_lists]
        list_leaves = list_of_lists_without_empty_strings
        list_lens = [len(i) for i in list_leaves]
        otu_train_cols = list()
        for i, j in zip(list_leaves, list_lens):
            if j == 1:
                updated = "k__" + i[0]
            elif j == 2:
                updated = "k__" + i[0] + ";" + "p__" + i[1]
            elif j == 3:
                updated = "k__" + i[0] + ";" + "p__" + i[1] + ";" + "c__" + \
                          i[2]
            elif j == 4:
                updated = "k__" + i[0] + ";" + "p__" + i[1] + ";" + "c__" + i[
                    2] + ";" + "o__" + i[3]
            elif j == 5:
                updated = "k__" + i[0] + ";" + "p__" + i[1] + ";" + "c__" + i[
                    2] + ";" + "o__" + i[3] + ";" + "f__" + i[4]
            elif j == 6:
                updated = "k__" + i[0] + ";" + "p__" + i[1] + ";" + "c__" + i[
                    2] + ";" + "o__" + i[3] + ";" + "f__" + i[4] + ";" + "g__" + \
                          i[5]
            elif j == 7:
                updated = "k__" + i[0] + ";" + "p__" + i[1] + ";" + "c__" + i[
                    2] + ";" + "o__" + i[3] + ";" + "f__" + i[4] + ";" + "g__" + i[
                              5] + ";" + "s__" + i[6]
            elif j == 8:
                updated = "k__" + i[0] + ";" + "p__" + i[1] + ";" + "c__" + i[
                    2] + ";" + "o__" + i[3] + ";" + "f__" + i[4] + ";" + "g__" + i[
                              5] + ";" + "s__" + i[6]
            elif j == 9:
                updated = "k__" + i[0] + ";" + "p__" + i[1] + ";" + "c__" + i[
                    2] + ";" + "o__" + i[3] + ";" + "f__" + i[4] + ";" + "g__" + i[
                              5] + ";" + "s__" + i[6]
            elif j == 10:
                updated = "k__" + i[0] + ";" + "p__" + i[1] + ";" + "c__" + i[
                    2] + ";" + "o__" + i[3] + ";" + "f__" + i[4] + ";" + "g__" + i[
                              5] + ";" + "s__" + i[6]
            elif j == 11:
                updated = "k__" + i[0] + ";" + "p__" + i[1] + ";" + "c__" + i[
                    2] + ";" + "o__" + i[3] + ";" + "f__" + i[4] + ";" + "g__" + i[
                              5] + ";" + "s__" + i[6]
            elif j == 12:
                updated = "k__" + i[0] + ";" + "p__" + i[1] + ";" + "c__" + i[
                    2] + ";" + "o__" + i[3] + ";" + "f__" + i[4] + ";" + "g__" + i[
                              5] + ";" + "s__" + i[6]
            otu_train_cols.append(updated)
        p.index = otu_train_cols
        p = p.groupby(p.index).first()
    elif model == "ALDEx Wilcoxon-C":
        p = pd.read_csv(f"ancom_bc2/data/{data}/aldex_real_holm.csv",index_col=0)[f"wi.eBH"]
        p.index = [i.replace("..", ";").replace(".", ";") for i in p.index]
        p.index = [i.replace("_p", ";p").replace("_c", ";c").replace("_o", ";o").replace("_f", ";f").replace("_g",
                                                                                                             ";g").replace(
            "_s", ";s") for i in p.index]
        if data == "MF":
            p.index = [i.split("t__")[0] for i in p.index]
        list_of_lists = [(i.replace("k__", "").replace("p__", "").replace("c__", "").replace("o__", "").replace("f__",
                                                                                                                "").replace(
            "g__", "").replace("s__", "")).split(";") for i in p.index]
        list_of_lists_without_empty_strings = [[item for item in sublist if item != ""] for sublist in list_of_lists]
        list_leaves = list_of_lists_without_empty_strings
        list_lens = [len(i) for i in list_leaves]
        otu_train_cols = list()
        for i, j in zip(list_leaves, list_lens):
            if j == 1:
                updated = "k__" + i[0]
            elif j == 2:
                updated = "k__" + i[0] + ";" + "p__" + i[1]
            elif j == 3:
                updated = "k__" + i[0] + ";" + "p__" + i[1] + ";" + "c__" + \
                          i[2]
            elif j == 4:
                updated = "k__" + i[0] + ";" + "p__" + i[1] + ";" + "c__" + i[
                    2] + ";" + "o__" + i[3]
            elif j == 5:
                updated = "k__" + i[0] + ";" + "p__" + i[1] + ";" + "c__" + i[
                    2] + ";" + "o__" + i[3] + ";" + "f__" + i[4]
            elif j == 6:
                updated = "k__" + i[0] + ";" + "p__" + i[1] + ";" + "c__" + i[
                    2] + ";" + "o__" + i[3] + ";" + "f__" + i[4] + ";" + "g__" + \
                          i[5]
            elif j == 7:
                updated = "k__" + i[0] + ";" + "p__" + i[1] + ";" + "c__" + i[
                    2] + ";" + "o__" + i[3] + ";" + "f__" + i[4] + ";" + "g__" + i[
                              5] + ";" + "s__" + i[6]
            elif j == 8:
                updated = "k__" + i[0] + ";" + "p__" + i[1] + ";" + "c__" + i[
                    2] + ";" + "o__" + i[3] + ";" + "f__" + i[4] + ";" + "g__" + i[
                              5] + ";" + "s__" + i[6]
            elif j == 9:
                updated = "k__" + i[0] + ";" + "p__" + i[1] + ";" + "c__" + i[
                    2] + ";" + "o__" + i[3] + ";" + "f__" + i[4] + ";" + "g__" + i[
                              5] + ";" + "s__" + i[6]
            elif j == 10:
                updated = "k__" + i[0] + ";" + "p__" + i[1] + ";" + "c__" + i[
                    2] + ";" + "o__" + i[3] + ";" + "f__" + i[4] + ";" + "g__" + i[
                              5] + ";" + "s__" + i[6]
            elif j == 11:
                updated = "k__" + i[0] + ";" + "p__" + i[1] + ";" + "c__" + i[
                    2] + ";" + "o__" + i[3] + ";" + "f__" + i[4] + ";" + "g__" + i[
                              5] + ";" + "s__" + i[6]
            elif j == 12:
                updated = "k__" + i[0] + ";" + "p__" + i[1] + ";" + "c__" + i[
                    2] + ";" + "o__" + i[3] + ";" + "f__" + i[4] + ";" + "g__" + i[
                              5] + ";" + "s__" + i[6]
            elif j == 13:
                updated = "k__" + i[0] + ";" + "p__" + i[1] + ";" + "c__" + i[
                    2] + ";" + "o__" + i[3] + ";" + "f__" + i[4] + ";" + "g__" + i[
                              5] + ";" + "s__" + i[6]
            elif j == 14:
                updated = "k__" + i[0] + ";" + "p__" + i[1] + ";" + "c__" + i[
                    2] + ";" + "o__" + i[3] + ";" + "f__" + i[4] + ";" + "g__" + i[
                              5] + ";" + "s__" + i[6]
            elif j == 15:
                updated = "k__" + i[0] + ";" + "p__" + i[1] + ";" + "c__" + i[
                    2] + ";" + "o__" + i[3] + ";" + "f__" + i[4] + ";" + "g__" + i[
                              5] + ";" + "s__" + i[6]

            otu_train_cols.append(updated)
        p.index = otu_train_cols
        p = p.groupby(p.index).first()
    elif model == "ANCOM-BC2":
        if DATA == "jacob":
            p = pd.read_csv(f"ancom_bc2/data/{data}/jacob_real_holm.csv", index_col=0)["p_Tag"]
            p.index = [i.replace("..", ";").replace(".", ";") for i in p.index]
            p.index = [i.replace("_p", ";p").replace("_c", ";c").replace("_o", ";o").replace("_f", ";f").replace("_g",
                                                                                                                 ";g").replace(
                "_s", ";s") for i in p.index]
            list_of_lists = [
                (i.replace("k__", "").replace("p__", "").replace("c__", "").replace("o__", "").replace("f__",
                                                                                                       "").replace(
                    "g__", "").replace("s__", "")).split(";") for i in p.index]
            list_of_lists_without_empty_strings = [[item for item in sublist if item != ""] for sublist in
                                                   list_of_lists]
            list_leaves = list_of_lists_without_empty_strings
            list_lens = [len(i) for i in list_leaves]
            otu_train_cols = list()
            for i, j in zip(list_leaves, list_lens):
                if j == 1:
                    updated = "k__" + i[0]
                elif j == 2:
                    updated = "k__" + i[0] + ";" + "p__" + i[1]
                elif j == 3:
                    updated = "k__" + i[0] + ";" + "p__" + i[1] + ";" + "c__" + \
                              i[2]
                elif j == 4:
                    updated = "k__" + i[0] + ";" + "p__" + i[1] + ";" + "c__" + i[
                        2] + ";" + "o__" + i[3]
                elif j == 5:
                    updated = "k__" + i[0] + ";" + "p__" + i[1] + ";" + "c__" + i[
                        2] + ";" + "o__" + i[3] + ";" + "f__" + i[4]
                elif j == 6:
                    updated = "k__" + i[0] + ";" + "p__" + i[1] + ";" + "c__" + i[
                        2] + ";" + "o__" + i[3] + ";" + "f__" + i[4] + ";" + "g__" + \
                              i[5]
                elif j == 7:
                    updated = "k__" + i[0] + ";" + "p__" + i[1] + ";" + "c__" + i[
                        2] + ";" + "o__" + i[3] + ";" + "f__" + i[4] + ";" + "g__" + i[
                                  5] + ";" + "s__" + i[6]
                elif j == 8:
                    updated = "k__" + i[0] + ";" + "p__" + i[1] + ";" + "c__" + i[
                        2] + ";" + "o__" + i[3] + ";" + "f__" + i[4] + ";" + "g__" + i[
                                  5] + ";" + "s__" + i[6]
                elif j == 9:
                    updated = "k__" + i[0] + ";" + "p__" + i[1] + ";" + "c__" + i[
                        2] + ";" + "o__" + i[3] + ";" + "f__" + i[4] + ";" + "g__" + i[
                                  5] + ";" + "s__" + i[6]
                elif j == 10:
                    updated = "k__" + i[0] + ";" + "p__" + i[1] + ";" + "c__" + i[
                        2] + ";" + "o__" + i[3] + ";" + "f__" + i[4] + ";" + "g__" + i[
                                  5] + ";" + "s__" + i[6]
                elif j == 11:
                    updated = "k__" + i[0] + ";" + "p__" + i[1] + ";" + "c__" + i[
                        2] + ";" + "o__" + i[3] + ";" + "f__" + i[4] + ";" + "g__" + i[
                                  5] + ";" + "s__" + i[6]
                elif j == 12:
                    updated = "k__" + i[0] + ";" + "p__" + i[1] + ";" + "c__" + i[
                        2] + ";" + "o__" + i[3] + ";" + "f__" + i[4] + ";" + "g__" + i[
                                  5] + ";" + "s__" + i[6]
                elif j == 13:
                    updated = "k__" + i[0] + ";" + "p__" + i[1] + ";" + "c__" + i[
                        2] + ";" + "o__" + i[3] + ";" + "f__" + i[4] + ";" + "g__" + i[
                                  5] + ";" + "s__" + i[6]
                elif j == 14:
                    updated = "k__" + i[0] + ";" + "p__" + i[1] + ";" + "c__" + i[
                        2] + ";" + "o__" + i[3] + ";" + "f__" + i[4] + ";" + "g__" + i[
                                  5] + ";" + "s__" + i[6]
                elif j == 15:
                    updated = "k__" + i[0] + ";" + "p__" + i[1] + ";" + "c__" + i[
                        2] + ";" + "o__" + i[3] + ";" + "f__" + i[4] + ";" + "g__" + i[
                                  5] + ";" + "s__" + i[6]
                otu_train_cols.append(updated)
            p.index = otu_train_cols
            p = p.groupby(p.index).first()
        else:
            p = pd.read_csv(f"ancom_bc2/data/{data}/real_holm.csv", index_col=0)["p_Tag"]
            p.index = [i.replace("..", ";").replace(".", ";") for i in p.index]

    elif model == "ANCOM-BC2-C":
        if DATA == "jacob":
            p = pd.read_csv(f"ancom_bc2/data/{data}/jacob_real_holm.csv", index_col=0)["q_Tag"]
            p.index = [i.replace("..", ";").replace(".", ";") for i in p.index]
            p.index = [i.replace("_p", ";p").replace("_c", ";c").replace("_o", ";o").replace("_f", ";f").replace("_g",
                                                                                                                 ";g").replace(
                "_s", ";s") for i in p.index]
            list_of_lists = [
                (i.replace("k__", "").replace("p__", "").replace("c__", "").replace("o__", "").replace("f__",
                                                                                                       "").replace(
                    "g__", "").replace("s__", "")).split(";") for i in p.index]
            list_of_lists_without_empty_strings = [[item for item in sublist if item != ""] for sublist in
                                                   list_of_lists]
            list_leaves = list_of_lists_without_empty_strings
            list_lens = [len(i) for i in list_leaves]
            otu_train_cols = list()
            for i, j in zip(list_leaves, list_lens):
                if j == 1:
                    updated = "k__" + i[0]
                elif j == 2:
                    updated = "k__" + i[0] + ";" + "p__" + i[1]
                elif j == 3:
                    updated = "k__" + i[0] + ";" + "p__" + i[1] + ";" + "c__" + \
                              i[2]
                elif j == 4:
                    updated = "k__" + i[0] + ";" + "p__" + i[1] + ";" + "c__" + i[
                        2] + ";" + "o__" + i[3]
                elif j == 5:
                    updated = "k__" + i[0] + ";" + "p__" + i[1] + ";" + "c__" + i[
                        2] + ";" + "o__" + i[3] + ";" + "f__" + i[4]
                elif j == 6:
                    updated = "k__" + i[0] + ";" + "p__" + i[1] + ";" + "c__" + i[
                        2] + ";" + "o__" + i[3] + ";" + "f__" + i[4] + ";" + "g__" + \
                              i[5]
                elif j == 7:
                    updated = "k__" + i[0] + ";" + "p__" + i[1] + ";" + "c__" + i[
                        2] + ";" + "o__" + i[3] + ";" + "f__" + i[4] + ";" + "g__" + i[
                                  5] + ";" + "s__" + i[6]
                elif j == 8:
                    updated = "k__" + i[0] + ";" + "p__" + i[1] + ";" + "c__" + i[
                        2] + ";" + "o__" + i[3] + ";" + "f__" + i[4] + ";" + "g__" + i[
                                  5] + ";" + "s__" + i[6]
                elif j == 9:
                    updated = "k__" + i[0] + ";" + "p__" + i[1] + ";" + "c__" + i[
                        2] + ";" + "o__" + i[3] + ";" + "f__" + i[4] + ";" + "g__" + i[
                                  5] + ";" + "s__" + i[6]
                elif j == 10:
                    updated = "k__" + i[0] + ";" + "p__" + i[1] + ";" + "c__" + i[
                        2] + ";" + "o__" + i[3] + ";" + "f__" + i[4] + ";" + "g__" + i[
                                  5] + ";" + "s__" + i[6]
                elif j == 11:
                    updated = "k__" + i[0] + ";" + "p__" + i[1] + ";" + "c__" + i[
                        2] + ";" + "o__" + i[3] + ";" + "f__" + i[4] + ";" + "g__" + i[
                                  5] + ";" + "s__" + i[6]
                elif j == 12:
                    updated = "k__" + i[0] + ";" + "p__" + i[1] + ";" + "c__" + i[
                        2] + ";" + "o__" + i[3] + ";" + "f__" + i[4] + ";" + "g__" + i[
                                  5] + ";" + "s__" + i[6]
                elif j == 13:
                    updated = "k__" + i[0] + ";" + "p__" + i[1] + ";" + "c__" + i[
                        2] + ";" + "o__" + i[3] + ";" + "f__" + i[4] + ";" + "g__" + i[
                                  5] + ";" + "s__" + i[6]
                elif j == 14:
                    updated = "k__" + i[0] + ";" + "p__" + i[1] + ";" + "c__" + i[
                        2] + ";" + "o__" + i[3] + ";" + "f__" + i[4] + ";" + "g__" + i[
                                  5] + ";" + "s__" + i[6]
                elif j == 15:
                    updated = "k__" + i[0] + ";" + "p__" + i[1] + ";" + "c__" + i[
                        2] + ";" + "o__" + i[3] + ";" + "f__" + i[4] + ";" + "g__" + i[
                                  5] + ";" + "s__" + i[6]
                otu_train_cols.append(updated)
            p.index = otu_train_cols
            p = p.groupby(p.index).first()

        else:
            p = pd.read_csv(f"ancom_bc2/data/{data}/real_holm.csv", index_col=0)["q_Tag"]
        p.index = [i.replace("..", ";").replace(".", ";") for i in p.index]

    elif model == "LEfSe":
        p = pd.read_csv(f"ancom_bc2/data/{data}/{data}_TRUE_lefse.csv", index_col=0)["P"]
        p = p[p.values != "-"]
        p.index = [i.replace("_p",";p").replace("_c",";c").replace("_o",";o").replace("_f",";f").replace("_g",";g").replace("_s",";s") for i in p.index]
        list_of_lists = [(i.replace("k__","").replace("p__","").replace("c__","").replace("o__","").replace("f__","").replace("g__","").replace("s__","")).split(";") for i in p.index]
        list_of_lists_without_empty_strings = [[item for item in sublist if item != ""] for sublist in list_of_lists]
        list_leaves = list_of_lists_without_empty_strings
        list_lens = [len(i) for i in list_leaves]
        otu_train_cols = list()
        for i, j in zip(list_leaves, list_lens):
            if j == 1:
                updated = "k__" + i[0]
            elif j == 2:
                updated = "k__" + i[0] + ";" + "p__" + i[1]
            elif j == 3:
                updated = "k__" + i[0] + ";" + "p__" + i[1] + ";" + "c__" + \
                              i[2]
            elif j == 4:
                updated = "k__" + i[0] + ";" + "p__" + i[1] + ";" + "c__" + i[
                        2] + ";" + "o__" + i[3]
            elif j == 5:
                updated = "k__" + i[0] + ";" + "p__" + i[1] + ";" + "c__" + i[
                        2] + ";" + "o__" + i[3] + ";" + "f__" + i[4]
            elif j == 6:
                updated = "k__" + i[0] + ";" + "p__" + i[1] + ";" + "c__" + i[
                        2] + ";" + "o__" + i[3] + ";" + "f__" + i[4] + ";" + "g__" + \
                              i[5]
            elif j == 7:
                updated = "k__" + i[0] + ";" + "p__" + i[1] + ";" + "c__" + i[
                        2] + ";" + "o__" + i[3] + ";" + "f__" + i[4] + ";" + "g__" + i[
                                  5] + ";" + "s__" + i[6]
            elif j == 8:
                updated = "k__" + i[0] + ";" + "p__" + i[1] + ";" + "c__" + i[
                        2] + ";" + "o__" + i[3] + ";" + "f__" + i[4] + ";" + "g__" + i[
                                  5] + ";" + "s__" + i[6]

            otu_train_cols.append(updated)
        p.index = otu_train_cols
        p = p.groupby(p.index).first()
        C=0
    elif model == "miMic":
        p = pd.read_csv(f"ancom_bc2/data/{data}/mimic_df_corrs.csv",index_col=0)["p"]
        if data == "Cirrhosis_Knight_Lab":
            p.index = [i[:-2] for i in p.index]
        list_leaves = list(p.index)
        list_lens = [len(i.split(";")) for i in list_leaves]
        otu_train_cols = list()
        if data == "Cirrhosis_Knight_Lab":
            for i, j in zip(list_leaves, list_lens):
                if j == 1:
                    updated = "k__" + i.split(";")[0]
                elif j == 2:
                    updated = "k__" + i.split(";")[0] + ";" + "p__" + i.split(";")[1]
                elif j == 3:
                    updated = "k__" + i.split(";")[0] + ";" + "p__" + i.split(";")[1] + ";" + "c__" + \
                              i.split(";")[2]
                elif j == 4:
                    updated = "k__" + i.split(";")[0] + ";" + "p__" + i.split(";")[1] + ";" + "c__" + i.split(";")[
                        2] + ";" + "o__" + i.split(";")[3]
                elif j == 5:
                    updated = "k__" + i.split(";")[0] + ";" + "p__" + i.split(";")[1] + ";" + "c__" + i.split(";")[
                        2] + ";" + "o__" + i.split(";")[3] + ";" + "f__" + i.split(";")[4]
                elif j == 6:
                    updated = "k__" + i.split(";")[0] + ";" + "p__" + i.split(";")[1] + ";" + "c__" + i.split(";")[
                        2] + ";" + "o__" + i.split(";")[3] + ";" + "f__" + i.split(";")[4] + ";" + "g__" + \
                              i.split(";")[5]
                elif j == 7:
                    updated = "k__" + i.split(";")[0] + ";" + "p__" + i.split(";")[1] + ";" + "c__" + i.split(";")[
                        2] + ";" + "o__" + i.split(";")[3] + ";" + "f__" + i.split(";")[4] + ";" + "g__" + i.split(";")[
                                  5] + ";" + "s__" + i.split(";")[6]
                otu_train_cols.append(updated)
        else:
            for i, j in zip(list_leaves, list_lens):
                if j == 1:
                    updated = "k__" + i.split(";")[0].split("_")[0]
                elif j == 2:
                    updated = "k__" + i.split(";")[0] + ";" + "p__" + i.split(";")[1].split("_")[0]
                elif j == 3:
                    updated = "k__" + i.split(";")[0] + ";" + "p__" + i.split(";")[1] + ";" + "c__" + \
                              i.split(";")[2].split("_")[0]
                elif j == 4:
                    updated = "k__" + i.split(";")[0] + ";" + "p__" + i.split(";")[1] + ";" + "c__" + i.split(";")[
                        2] + ";" + "o__" + i.split(";")[3].split("_")[0]
                elif j == 5:
                    updated = "k__" + i.split(";")[0] + ";" + "p__" + i.split(";")[1] + ";" + "c__" + i.split(";")[
                        2] + ";" + "o__" + i.split(";")[3] + ";" + "f__" + i.split(";")[4].split("_")[0]
                elif j == 6:
                    updated = "k__" + i.split(";")[0] + ";" + "p__" + i.split(";")[1] + ";" + "c__" + i.split(";")[
                        2] + ";" + "o__" + i.split(";")[3] + ";" + "f__" + i.split(";")[4] + ";" + "g__" + \
                              i.split(";")[5].split("_")[0]
                elif j == 7:
                    updated = "k__" + i.split(";")[0] + ";" + "p__" + i.split(";")[1] + ";" + "c__" + i.split(";")[
                        2] + ";" + "o__" + i.split(";")[3] + ";" + "f__" + i.split(";")[4] + ";" + "g__" + i.split(";")[
                                  5] + ";" + "s__" + i.split(";")[6].split("_")[0]
                otu_train_cols.append(updated)
        p.index = otu_train_cols

    elif model == "miMic relative":
        p = pd.read_csv(f"ancom_bc2/data/{data}/mimic_df_corrs_relative.csv",index_col=0)["p"]
        list_leaves = list(p.index)
        list_lens = [len(i.split(";")) for i in list_leaves]
        otu_train_cols = list()
        if data == "Cirrhosis_Knight_Lab":
            for i, j in zip(list_leaves, list_lens):
                if j == 1:
                    updated = "k__" + i.split(";")[0]
                elif j == 2:
                    updated = "k__" + i.split(";")[0] + ";" + "p__" + i.split(";")[1]
                elif j == 3:
                    updated = "k__" + i.split(";")[0] + ";" + "p__" + i.split(";")[1] + ";" + "c__" + \
                              i.split(";")[2]
                elif j == 4:
                    updated = "k__" + i.split(";")[0] + ";" + "p__" + i.split(";")[1] + ";" + "c__" + i.split(";")[
                        2] + ";" + "o__" + i.split(";")[3]
                elif j == 5:
                    updated = "k__" + i.split(";")[0] + ";" + "p__" + i.split(";")[1] + ";" + "c__" + i.split(";")[
                        2] + ";" + "o__" + i.split(";")[3] + ";" + "f__" + i.split(";")[4]
                elif j == 6:
                    updated = "k__" + i.split(";")[0] + ";" + "p__" + i.split(";")[1] + ";" + "c__" + i.split(";")[
                        2] + ";" + "o__" + i.split(";")[3] + ";" + "f__" + i.split(";")[4] + ";" + "g__" + \
                              i.split(";")[5]
                elif j == 7:
                    updated = "k__" + i.split(";")[0] + ";" + "p__" + i.split(";")[1] + ";" + "c__" + i.split(";")[
                        2] + ";" + "o__" + i.split(";")[3] + ";" + "f__" + i.split(";")[4] + ";" + "g__" + i.split(";")[
                                  5] + ";" + "s__" + i.split(";")[6]
                otu_train_cols.append(updated)
        else:
            for i, j in zip(list_leaves, list_lens):
                if j == 1:
                    updated = "k__" + i.split(";")[0].split("_")[0]
                elif j == 2:
                    updated = "k__" + i.split(";")[0] + ";" + "p__" + i.split(";")[1].split("_")[0]
                elif j == 3:
                    updated = "k__" + i.split(";")[0] + ";" + "p__" + i.split(";")[1] + ";" + "c__" + \
                              i.split(";")[2].split("_")[0]
                elif j == 4:
                    updated = "k__" + i.split(";")[0] + ";" + "p__" + i.split(";")[1] + ";" + "c__" + i.split(";")[
                        2] + ";" + "o__" + i.split(";")[3].split("_")[0]
                elif j == 5:
                    updated = "k__" + i.split(";")[0] + ";" + "p__" + i.split(";")[1] + ";" + "c__" + i.split(";")[
                        2] + ";" + "o__" + i.split(";")[3] + ";" + "f__" + i.split(";")[4].split("_")[0]
                elif j == 6:
                    updated = "k__" + i.split(";")[0] + ";" + "p__" + i.split(";")[1] + ";" + "c__" + i.split(";")[
                        2] + ";" + "o__" + i.split(";")[3] + ";" + "f__" + i.split(";")[4] + ";" + "g__" + \
                              i.split(";")[5].split("_")[0]
                elif j == 7:
                    updated = "k__" + i.split(";")[0] + ";" + "p__" + i.split(";")[1] + ";" + "c__" + i.split(";")[
                        2] + ";" + "o__" + i.split(";")[3] + ";" + "f__" + i.split(";")[4] + ";" + "g__" + i.split(";")[
                                  5] + ";" + "s__" + i.split(";")[6].split("_")[0]
                otu_train_cols.append(updated)
        p.index = otu_train_cols

        p = p.groupby(p.index).first()
    return p

def adjust_p_vals(p,MODEL):
    if isinstance(p, pd.Series):
        p = p.to_frame(MODEL)
    elif isinstance(p, pd.DataFrame):
        p.columns = [MODEL]
    return p

if __name__ == '__main__':
    MERGE = True
    mpl.rc('font', family='Times New Roman')
    LIST_COLORS = ["mediumorchid","hotpink","lightskyblue", "cornflowerblue", "salmon", "tomato", "mediumaquamarine", "seagreen", "grey",
                   "dimgrey", "lightslategray", "slategray", "green", "darkgreen", 'gold']

    LIST_MODELS = ["miMic relative","miMic","LINDA", "LINDA-C", "DeSeq", "DeSeq-C", "ANCOM", "ANCOM-C", "ALDEx Welch", "ALDEx Welch-C",
                   "ALDEx Wilcoxon", "ALDEx Wilcoxon-C", "ANCOM-BC2", "ANCOM-BC2-C", "LEfSe"]
    LIST_DATASETS_16S = ["Cirrhosis_Knight_Lab", "ERP020401", "ERP021216", "IBD", "he", "jacob", "ok_30", "ok_94",
                         "PRJNA353587", "PRJNA419097", "WB","MF"]

    LIST_DATASETS_WGS = ["Cirrhosis", "data1_metab", "data2_metab", "data3_metab", "data5_metab", "data6_metab", "T2D",
                         "obesity"]
    list_all_dfs = []

    for DATA in LIST_DATASETS_16S:
        print(DATA)
        significant_per_data = pd.DataFrame(columns=LIST_MODELS)
        list_total_asvs = []
        for MODEL in LIST_MODELS:
            print(MODEL)
            if MODEL == "LINDA" and DATA == "IBD":
                continue
            elif MODEL == "LINDA-C" and DATA == "IBD":
                continue
            elif MODEL == "LINDA-C" and DATA == "MF":
                continue
            elif MODEL == "LINDA" and DATA == "MF":
                continue
            elif MODEL == "LEfSe" and DATA == "MF":
                continue
            elif MODEL == "miMic" and DATA == "MF":
                continue
            elif MODEL == "miMic relative" and DATA == "MF":
                continue
            elif MODEL == "miMic relative" and DATA == "PRJNA353587":
                continue
            elif MODEL == "miMic relative" and DATA == "PRJNA419097":
                continue
            elif MODEL == "miMic relative" and DATA == "ERP020401":
                continue



            else:
                p = load_p_vals(MODEL,DATA)
                if MODEL == "ANCOM" or MODEL == "ANCOM-C" or MODEL == "LEfSe":
                    p = p.to_frame(MODEL)
                    significant = p


                else:
                    p = adjust_p_vals(p,MODEL)
                    significant = p[p.values < 0.05]
                list_total_asvs.append(significant)
        all_models = pd.concat(list_total_asvs,axis=1)
        all_models = all_models.fillna(2.0)
        all_models = all_models!=2.0
        all_models["Total"] = all_models.sum(axis=1)

        all_models.index = [f"{DATA}_{i}" for i in range(len(all_models.index))]
        list_all_dfs.append(all_models)
    # merge to 1 table for all datasets
    merged = pd.concat(list_all_dfs)

    # plot
    fig, axs = plt.subplots(nrows=5, ncols=3, figsize=(15, 15), sharex='col', sharey='row')
    list_mm = []
    for (MODEL, COLOR) in zip(LIST_MODELS, LIST_COLORS):
        if MERGE:
            if MODEL == "miMic relative":
                mm = build_combined_df(merged, LIST_MODELS)
                mm = mm[[MODEL, "b-LINDA", "b-ANCOM", "b-DeSeq", "b-ANCOM-BC2", "b-ALDEx", "LEfSe"]]
            elif MODEL == "miMic":
                mm = build_combined_df(merged, LIST_MODELS)
                mm = mm[[MODEL, "b-LINDA", "b-ANCOM", "b-DeSeq", "b-ANCOM-BC2", "b-ALDEx", "LEfSe"]]
            elif MODEL == "DeSeq":
                mm = build_combined_df(merged, LIST_MODELS)
                mm = mm[[MODEL, "b-LINDA", "b-ANCOM", "b-miMic", "b-ANCOM-BC2", "b-ALDEx", "LEfSe"]]
            elif MODEL == "DeSeq-C":
                mm = build_combined_df(merged, LIST_MODELS)
                mm = mm[[MODEL, "b-LINDA", "b-ANCOM", "b-miMic", "b-ANCOM-BC2", "b-ALDEx", "LEfSe"]]
            elif MODEL == "ANCOM":
                mm = build_combined_df(merged, LIST_MODELS)
                mm = mm[[MODEL, "b-LINDA", "b-DeSeq", "b-miMic", "b-ANCOM-BC2", "b-ALDEx", "LEfSe"]]
            elif MODEL == "ANCOM-C":
                mm = build_combined_df(merged, LIST_MODELS)
                mm = mm[[MODEL, "b-LINDA", "b-DeSeq", "b-miMic", "b-ANCOM-BC2", "b-ALDEx", "LEfSe"]]
            elif MODEL == "ANCOM-BC2":
                mm = build_combined_df(merged, LIST_MODELS)
                mm = mm[[MODEL, "b-LINDA", "b-DeSeq", "b-miMic", "b-ANCOM", "b-ALDEx", "LEfSe"]]
            elif MODEL == "ANCOM-BC2-C":
                mm = build_combined_df(merged, LIST_MODELS)
                mm = mm[[MODEL, "b-LINDA", "b-DeSeq", "b-miMic", "b-ANCOM", "b-ALDEx", "LEfSe"]]
            elif MODEL == "LEfSe":
                mm = build_combined_df(merged, LIST_MODELS)
                mm = mm[[MODEL, "b-LINDA", "b-DeSeq", "b-miMic", "b-ANCOM", "b-ALDEx", "ANCOM-BC2"]]
            elif MODEL == "ALDEx Welch":
                mm = build_combined_df(merged, LIST_MODELS)
                mm = mm[[MODEL, "b-LINDA", "b-DeSeq", "b-miMic", "b-ANCOM", "LEfSe", "ANCOM-BC2"]]
            elif MODEL == "ALDEx Welch-C":
                mm = build_combined_df(merged, LIST_MODELS)
                mm = mm[[MODEL, "b-LINDA", "b-DeSeq", "b-miMic", "b-ANCOM", "LEfSe", "ANCOM-BC2"]]
            elif MODEL == "ALDEx Wilcoxon":
                mm = build_combined_df(merged, LIST_MODELS)
                mm = mm[[MODEL, "b-LINDA", "b-DeSeq", "b-miMic", "b-ANCOM", "LEfSe", "ANCOM-BC2"]]
            elif MODEL == "ALDEx Wilcoxon-C":
                mm = build_combined_df(merged, LIST_MODELS)
                mm = mm[[MODEL, "b-LINDA", "b-DeSeq", "b-miMic", "b-ANCOM", "LEfSe", "ANCOM-BC2"]]
            elif MODEL == "LINDA":
                mm = build_combined_df(merged, LIST_MODELS)
                mm = mm[[MODEL, "b-ALDEx", "b-DeSeq", "b-miMic", "b-ANCOM", "LEfSe", "ANCOM-BC2"]]
            elif MODEL == "LINDA-C":
                mm = build_combined_df(merged, LIST_MODELS)
                mm = mm[[MODEL, "b-ALDEx", "b-DeSeq", "b-miMic", "b-ANCOM", "LEfSe", "ANCOM-BC2"]]
            mm = mm.fillna(False)
            mm["Total"] = mm.sum(axis=1)
        counts_df = plot_consistency_model(mm,MODEL,COLOR,MERGE)
        if MODEL == "miMic":
            mimic_counts_df = counts_df
        if MODEL not in ["miMic relative","miMic"]:
            list_mm.append(counts_df)

    merged = pd.concat(list_mm)
    avg_to_plot = merged.groupby(merged.index).mean()
    plot_paired_consistency_model(avg_to_plot,mimic_counts_df,"white","hotpink")




    c=0





import tqdm
from samba import micro2matrix
from scipy.stats import spearmanr, ttest_ind
import statsmodels.stats.multitest as smt
import pandas as pd
from scipy.stats import spearmanr
import os
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns


def treat_binary(tag_b):
    for col in tag_b.columns:
        if col != "gender":
            tag_b[col][tag_b[col] == "n"] = 0.0
            tag_b[col][tag_b[col] == "y"] = 1.0
            tag_b[col][tag_b[col] == "none"] = 2.0
        else:
            tag_b[col][tag_b[col] == "m"] = 0.0
            tag_b[col][tag_b[col] == "f"] = 1.0
    return tag_b


def treat_continuous(tag_c):
    for col in tag_c.columns:
        print(col)
        med = tag_c[col][tag_c[col] != "none"].median()
        tag_c[col] = tag_c[col].replace('none', med)
        tag_c[col] = tag_c[col].astype(float)
        mean = tag_c[col].mean()
        std = tag_c[col].std()
        tag_c[col] = (tag_c[col] - mean) / std
    return tag_c


def check_corr(df, meta, name):
    num_s_df = pd.DataFrame(index=meta.columns, columns=["number"])
    for col in meta.columns:
        corrs_df = pd.DataFrame(index=df.columns, columns=["scc", "p"])
        for bact in df.columns:
            scc, p = spearmanr(df[bact], meta[col])
            corrs_df["scc"][bact] = scc
            corrs_df["p"][bact] = p
        corrs_df.to_csv(f"metric_results/new/ggmp_{name}_shuffle/{col}.csv")

        n_significant = (corrs_df["p"] < 0.05).sum()
        num_s_df["number"][col] = n_significant
    num_s_df.to_csv(f"metric_results/new/ggmp_{name}_shuffle/number_significant.csv")


def load_img(folder_path, tag):
    arrays = []
    names = []
    for file in os.listdir(folder_path):
        if file.endswith(".npy"):
            if file == "bact_names.npy":
                continue
            file_path = os.path.join(folder_path, file)
            if file_path.split("\\")[-1].replace(".npy", "") in [str(i) for i in tag.index]:
                arrays.append(np.load(file_path, allow_pickle=True, mmap_mode='r'))
                names.append(file_path.split("\\")[-1].replace(".npy", ""))
    # final_array = np.concatenate(arrays)
    final_array = np.stack(arrays, axis=0)
    return final_array, names


def get_row_and_col(bact_df, taxon):
    row_index, col_index = bact_df.index[bact_df.eq(taxon).any(axis=1)][0], bact_df.columns[bact_df.eq(taxon).any()]
    return col_index, row_index


def calc_unique_corr(bact_df, taxon, imgs, tag,eval="corr"):
    col_index, row_index = get_row_and_col(bact_df, taxon)
    first_col_index = list(col_index)[0]
    result = [imgs[i, row_index, first_col_index] for i in range(imgs.shape[0])]
    # calc corr
    if eval == "corr":
        scc, p = spearmanr(result, tag)
    else:
        zero_indexes = tag.index[tag == 0]
        one_indexes = tag.index[tag == 1]
        otu0 = [result[i] for i in zero_indexes.tolist()]
        otu1 = [result[i] for i in one_indexes.tolist()]
        scc,p =ttest_ind(otu0,otu1)
    return col_index, row_index, scc, p


def calculate_all_imgs_tag_corr(mat, names, tag, start_i, algo="weighted",eval="corr",sis=None):
    img_arrays = mat
    bact_names = names  # np.load(f'{folder}/bact_names.npy', allow_pickle=True)
    bact_names_df = pd.DataFrame(bact_names)

    dict_corrs = dict()
    dict_ps = dict()
    all_ps = dict()
    different_tax_in_level = list(set(bact_names_df.iloc[start_i]))  # tax 1

    def binary_rec_by_pval(different_tax_in_level, eval="corr", sis=None):
        for tax in different_tax_in_level:
            # Stop condition
            if tax == '' or tax == "0.0" or tax is None:
                # Leaf in the middle of the tree
                return "leaf"

            col_index, row_index, scc, p = calc_unique_corr(bact_names_df, tax, img_arrays, tag, eval)
            all_ps[tax] = p
            if p >= 0.05:
                # Not significant - stop
                continue
            if (row_index + 1) >= bact_names_df.shape[0]:
                # Leaf in the end of the tree
                dict_corrs[tax] = scc
                dict_ps[tax] = p
                continue

            all_sons = set(bact_names_df[col_index].loc[row_index + 1])
            ret = binary_rec_by_pval(all_sons, eval,sis)
            if ret == "leaf":
                # This is the leaf:
                dict_corrs[tax] = scc
                dict_ps[tax] = p

            if sis == "bonferroni" and len(all_sons) > 1 and len(all_sons.intersection(dict_ps.keys())) > 0:
                sons_pv = {k: all_ps[k] for k in all_sons}
                min_son = min(sons_pv, key=sons_pv.get)
                del sons_pv[min_son]

                rejected_r, corrected_p_values_r, _, _ = smt.multipletests(list(sons_pv.values()),
                                                                           method="bonferroni")
                for e, son in enumerate(sons_pv):
                    if corrected_p_values_r[e] >= 0.05:
                        for bact in [k for k in dict_ps.keys() if son in k]:
                            del dict_ps[bact]
                    else:
                        dict_ps[son] = corrected_p_values_r[e]

    def rec_weighted_update(different_tax_in_level, val_pred=0, eval="corr"):
        for tax in different_tax_in_level:
            # Stop condition
            if tax == '':
                # Leaf in the middle of the tree
                return "leaf"

            col_index, row_index = get_row_and_col(bact_names_df, tax)

            if type(val_pred) != int:
                img_arrays[:, row_index, col_index] += val_pred.reshape(-1, 1)

            if (row_index + 1) >= bact_names_df.shape[0]:
                # Leaf in the end of the tree
                continue

            p = calc_unique_corr(bact_names_df, tax, img_arrays, tag, eval)[-1]
            all_sons = set(bact_names_df[col_index].loc[row_index + 1])

            if p >= 0.05:
                ret = rec_weighted_update(all_sons, val_pred, eval)
            else:
                ret = rec_weighted_update(all_sons, img_arrays[:, row_index, col_index[0]] / len(all_sons), eval)

            if ret == "leaf":
                # This is the leaf:
                if type(val_pred) != int:
                    img_arrays[:, row_index, col_index] += val_pred.reshape(-1, 1)

    def rec_all_leafs(different_tax_in_level, eval):
        for tax in different_tax_in_level:
            # Stop condition
            if tax == '':
                # Leaf in the middle of the tree
                return "leaf"

            col_index, row_index, scc, p = calc_unique_corr(bact_names_df, tax, img_arrays, tag, eval)
            if (row_index + 1) >= bact_names_df.shape[0]:
                # Leaf in the end of the tree
                dict_corrs[tax] = scc
                dict_ps[tax] = p
                continue

            all_sons = set(bact_names_df[col_index].loc[row_index + 1])
            ret = rec_all_leafs(all_sons, eval)
            if ret == "leaf":
                # This is the leaf:
                dict_corrs[tax] = scc
                dict_ps[tax] = p

    if algo == "binary":
        binary_rec_by_pval(different_tax_in_level, eval,sis)
    elif algo == "weighted":
        rec_weighted_update(different_tax_in_level, eval)
        rec_all_leafs(different_tax_in_level, eval)
    else:
        raise Exception()

    series_corrs = pd.Series(dict_corrs).to_frame("scc")
    series_ps = pd.Series(dict_ps).to_frame("p")
    df_corrs = pd.concat([series_corrs, series_ps], axis=1)
    df_corrs.index = [i.split("_")[0] for i in df_corrs.index]
    df_corrs = df_corrs.groupby(df_corrs.index).mean()
    return df_corrs


if __name__ == '__main__':
    # CHOOSE 1 OF {"img","1d"}
    KIND = "img"
    # CHOOSE 1 OF {"ttest","corr"}
    EVAL = "ttest" #"corr"
    mean1 = 0  # Mean of the normal distribution
    std1 = 1
    RUNS = 100
    list_sizes = [20, 40, 80, 160, 320, 640, 1280]
    list_mu = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
    TP_TABLE = pd.DataFrame(index=list_sizes, columns=list_mu)
    FP_TABLE = pd.DataFrame(index=list_sizes, columns=list_mu)
    tqdm_ = tqdm.tqdm(total=len(list_sizes) * len(list_mu) * RUNS)
    n = 0
    for size in list_sizes:  # Size of the vector
        for mean2 in list_mu:
            number_tp = 0
            number_fp = 0
            for run in range(RUNS):
                tqdm_.update()
                half_size = size / 2
                # Generate a vector of size 100 with values distributed normally
                b1 = np.random.normal(mean1, std1, size)
                b2 = np.random.normal(mean1, std1, size)
                b3_a = np.random.normal(mean1, std1, int(half_size))
                b3_b = np.random.normal(mean2, std1, int(half_size))
                t_a = [0.0 for i in range(int(half_size))]
                t_b = [1.0 for i in range(int(half_size))]
                t = t_a + t_b
                b3 = np.concatenate((b3_a, b3_b))

                otu_data = pd.DataFrame(index=range(size),
                                        columns=["k__Bacteria;p__Aa", "k__Bacteria;p__Bb", "k__Bacteria;p__Cc"])
                otu_data["k__Bacteria;p__Aa"] = b1
                otu_data["k__Bacteria;p__Bb"] = b2
                otu_data["k__Bacteria;p__Cc"] = b3
                tag_data = pd.DataFrame(index=otu_data.index, columns=["Tag"])
                tag_data["Tag"] = t
                if KIND != "img":
                    # calc raw corr
                    corrs_df = pd.DataFrame(index=otu_data.columns, columns=["scc", "p"])
                    for col in otu_data.columns:
                        if EVAL == "corr":
                            scc, p = spearmanr(otu_data[col], tag_data["Tag"])
                            corrs_df["scc"][col] = scc
                            corrs_df["p"][col] = p
                        else:
                            tag1 = tag_data[tag_data["Tag"] ==1.0]
                            tag0 = tag_data[tag_data["Tag"] ==0.0]
                            otu1 = otu_data.loc[tag1.index]
                            otu0 = otu_data.loc[tag0.index]
                            scc,p = ttest_ind(otu0[col],otu1[col])
                            corrs_df["scc"][col] = scc
                            corrs_df["p"][col] = p
                    significant = corrs_df[corrs_df["p"] < 0.05]
                    if 'k__Bacteria;p__Cc' in significant.index:
                        number_tp = number_tp + 1
                    if 'k__Bacteria;p__Aa' in significant.index or 'k__Bacteria;p__Bb' in significant.index:
                        number_fp = number_fp + 1
                else:

                    mat, names = micro2matrix(otu_data, folder=f"imgs_{size}_{mean2}", save=False)
                    TASK = "binary"
                    #
                    df_corrs = df_corrs = calculate_all_imgs_tag_corr(mat, names, tag_data["Tag"], 1, algo=TASK,eval=EVAL,sis="bonferroni")
                    if 'Bacteria;Cc' in df_corrs.index:
                        number_tp = number_tp + 1
                    if 'Bacteria;Aa' in df_corrs.index or 'Bacteria;Bb' in df_corrs.index:
                        number_fp = number_fp + 1
            FP_TABLE[mean2][size] = number_fp
            TP_TABLE[mean2][size] = number_tp
        if KIND == "img":
            FP_TABLE.to_csv(f"RESULTS/OURS/fp_table_0_0_mu_sis.csv")
            TP_TABLE.to_csv(f"RESULTS/OURS/tp_table_0_0_mu_sis.csv")
        else:
            FP_TABLE.to_csv(f"RESULTS/1D/fp_table_t.csv")
            TP_TABLE.to_csv(f"RESULTS/1D/tp_table_t.csv")

    plt.figure(figsize=(8, 6))  # Adjust the figure size as needed
    sns.heatmap(TP_TABLE.astype(float), annot=True, cmap='rocket_r')
    plt.title("TP")
    plt.show()
    plt.clf()

    plt.figure(figsize=(8, 6))  # Adjust the figure size as needed
    sns.heatmap(FP_TABLE.astype(float), annot=True, cmap='rocket_r')
    plt.title("FP")
    plt.show()
    plt.clf()

    c = 0

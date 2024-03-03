import pandas as pd
import numpy as np
import itertools
from itertools import chain
from scipy.stats import spearmanr
import matplotlib as mpl
import matplotlib.pyplot as plt

def load_data(folder, p, shuffle=False):
    if folder == "data1_metab":
        otu = pd.read_csv(f"Data_for_git_process/{folder}/tax7_relative_mean.csv", index_col=0)
        tags = pd.read_csv(f"Data_for_git_process/{folder}/tag_for_learning.csv", index_col=0)["Tag"]
        folder_ = f"2D_otus_dendogram_ordered/{folder}"
    elif folder == "data2_metab":
        otu = pd.read_csv(f"Data_for_git_process/{folder}/tax7_relative.csv", index_col=0)
        otu.index = [str(i) for i in otu.index]
        tags = pd.read_csv(f"Data_for_git_process/{folder}/tag.csv", index_col=0)
        tags.index = [str(i) for i in tags.index]
        folder_ = f"2D_otus_dendogram_ordered/{folder}"

    elif folder == "data6_metab":
        otu = pd.read_csv(f"Data_for_git_process/{folder}/tax7_relative_mean.csv", index_col=0)
        tags = pd.read_csv(f"Data_for_git_process/{folder}/tag.csv", index_col=0)
        folder_ = f"2D_otus_dendogram_ordered/{folder}"

    elif folder == "data5_metab":
        otu = pd.read_csv(f"Data_for_git_process/{folder}/tax7_relative_mean.csv", index_col=0)
        tags = pd.read_csv(f"Data_for_git_process/{folder}/tag.csv", index_col=0)
        tags0 = tags[tags["Study.Group"] == 0.0]
        tags1 = tags[tags["Study.Group"] == 1.0]
        folder_ = f"2D_otus_dendogram_ordered/{folder}"

    elif folder == "data3_metab":
        otu = pd.read_csv(f"Data_for_git_process/{folder}/tax7_relative_mean.csv", index_col=0)
        tags = pd.read_csv(f"Data_for_git_process/{folder}/tag.csv", index_col=0)
        tags0 = tags[tags["Study.Group"] == 0.0]
        tags1 = tags[tags["Study.Group"] == 1.0]
        folder_ = f"2D_otus_dendogram_ordered/{folder}"

    elif folder == "Cirrhosis_Knight_Lab":
        otu = pd.read_csv(f"Data/{folder}/mean_tax_7_relative_rare_bact_5_without_viruses.csv", index_col=0)
        tags = pd.read_csv(f"Data/{folder}/tag.csv", index_col=0)
        tags0 = tags[tags["Var"] == 0.0]
        tags1 = tags[tags["Var"] == 1.0]
        folder_ = f"2D_otus_dendogram_ordered/Cirrhosis_new"

    elif folder == "White_vs_black_vagina":
        otu = pd.read_csv(f"Data/{folder}/sum_relative.csv", index_col=0)
        tags = pd.read_csv(f"Data/{folder}/tag.csv", index_col=0)
        tags0 = tags[tags["Var"] == 0.0]
        tags1 = tags[tags["Var"] == 1.0]
        folder_ = f"2D_otus_dendogram_ordered/{folder}"

    elif folder == "Male_vs_female":
        otu = pd.read_csv(f"Data/{folder}/sub_pca_tax_7_log_rare_bact_5_without_viruses.csv", index_col=0)
        tags = pd.read_csv(f"Data/{folder}/tag.csv", index_col=0)
        tags0 = tags[tags["Var"] == 0.0]
        tags1 = tags[tags["Var"] == 1.0]
        folder_ = f"2D_otus_dendogram_ordered/MF_new"

    elif folder == "Nugent":
        otu = pd.read_csv(f"Data/White_vs_black_vagina/sum_relative.csv", index_col=0)
        tags = pd.read_csv(f"Data/Nugent_vagina/tag1.csv", index_col=0)
        tags0 = tags[tags["Tag"] == 0.0]
        tags1 = tags[tags["Tag"] == 1.0]

    elif folder == "Allergy/New_allergy_milk":
        otu = pd.read_csv(f"Data/Allergy/sum_relative.csv", index_col=0)
        tags = pd.read_csv(f"Data/{folder}/tag.csv", index_col=0)
        tags0 = tags[tags["Tag"] == 0.0]
        tags1 = tags[tags["Tag"] == 1.0]


    elif folder == "Allergy/New_allergy_nuts":
        otu = pd.read_csv(f"Data/Allergy/sum_relative.csv", index_col=0)
        tags = pd.read_csv(f"Data/{folder}/tag.csv", index_col=0)
        tags0 = tags[tags["AllergyType"] == 0.0]
        tags1 = tags[tags["AllergyType"] == 1.0]

    elif folder == "Allergy/New_allergy_Peanuts":
        otu = pd.read_csv(f"Data/Allergy/sum_relative.csv", index_col=0)
        tags = pd.read_csv(f"Data/{folder}/tag.csv", index_col=0)
        tags0 = tags[tags["AllergyType"] == 0.0]
        tags1 = tags[tags["AllergyType"] == 1.0]

    elif folder == "GDM_OLD":
        otu = pd.read_csv(f"Data/{folder}/otu.csv", index_col=0)
        tags = pd.read_csv(f"Data/{folder}/tag.csv", index_col=0)
        tags0 = tags[tags["Status"] == 0.0]
        tags1 = tags[tags["Status"] == 1.0]

    elif folder == "jacob":
        otu = pd.read_csv(f"Datasets_for_thesis_2_preprocess/efrat_data/after/{folder}_otu_relative_mean_tax7.csv",
                          index_col=0)
        #otu = pd.read_csv("Datasets_for_thesis_2_preprocess/efrat_data/after/jacob_otu_sub_pca_log_zscore_row_tax7.csv",index_col=0)
        tags = pd.read_csv(f"Data_for_git_process/{folder}/tag.csv", index_col=0)
        folder_ = f"2D_otus_dendogram_ordered/{folder}"

    elif folder == "kim":
        otu = pd.read_csv(f"Datasets_for_thesis_2_preprocess/efrat_data/after/{folder}_otu_relative_mean_tax7.csv",
                          index_col=0)
        tags = pd.read_csv(f"Data_for_git_process/{folder}/tag.csv", index_col=0)
        folder_ = f"2D_otus_dendogram_ordered/{folder}"

    elif folder == "he":
        otu = pd.read_csv(f"Datasets_for_thesis_2_preprocess/efrat_data/after/he_tax7_relative_mean.csv", index_col=0)
        #otu = pd.read_csv("Datasets_for_thesis_2_preprocess/efrat_data/after/he_tax7_subpca_zscore_row.csv",index_col=0)
        tags = pd.read_csv(f"Data_for_git_process/{folder}/tag.csv", index_col=0)
        tags.index = ["12021." + i.replace("-", ".") for i in tags.index]
        folder_ = f"2D_otus_dendogram_ordered/{folder}"

    elif folder == "T2D":
        otu = pd.read_csv(
            f"Data/PopPhy_data/{folder}/T2D_otu_mean_relative_rare_bact_5_tax_7_only_after_git_preprocess.csv",
            index_col=0)
        otu.index = [str(i) for i in otu.index]
        tags = pd.read_csv(f"Data/PopPhy_data/{folder}/tag.csv", index_col=0)
        tags.index = [str(i) for i in tags.index]
        folder_ = f"2D_otus_dendogram_ordered/{folder}"

    elif folder == "Obesity":
        otu = pd.read_csv(f"Data/PopPhy_data/{folder}/log_subpca_tax_7.csv", index_col=0)
        otu.index = [str(i) for i in otu.index]
        tags = pd.read_csv(f"Data/PopPhy_data/{folder}/tag.csv", index_col=0)["Tag"]
        tags.index = [str(i) for i in tags.index]
        folder_ = f"2D_otus_dendogram_ordered/{folder}"

    elif folder == "Cirrhosis":
        otu = pd.read_csv(f"Data/PopPhy_data/{folder}/log_subpca_tax_7.csv", index_col=0)
        otu.index = [str(i) for i in otu.index]
        tags = pd.read_csv(f"Data/PopPhy_data/{folder}/tag.csv", index_col=0)["Tag"]
        tags.index = [str(i) for i in tags.index]
        folder_ = f"2D_otus_dendogram_ordered/Cirrhosis_pop"

    elif folder == "IBD":
        #otu = pd.read_csv(f"Data/{folder}/only_otus_sum_relative_tax_7_rare_5.csv", index_col=0)
        otu = pd.read_csv(f"Data/IBD/only_otus_sub_pca_log_tax_7_rare_5.csv", index_col=0)
        otu.index = [str(i) for i in otu.index]
        tags = pd.read_csv(f"Data/{folder}/tag_IBD_VS_ALL.csv", index_col=0)["Tag"]
        tags.index = [str(i) for i in tags.index]
        folder_ = f"2D_otus_dendogram_ordered/{folder}"

    elif folder == "ok_108-SIBO":
        otu = pd.read_csv(f"Data/ok_108/tax7_log_subpca.csv", index_col=0)
        otu.index = [str(i) for i in otu.index]
        tags = pd.read_csv(f"Raw_data/ok_108/tag_SIBO.csv", index_col=0)["Tag"]
        tags.index = [str(i) for i in tags.index]

    elif folder == "ok_108-MF":
        otu = pd.read_csv(f"Data/ok_108/tax7_log_subpca.csv", index_col=0)
        otu.index = [str(i) for i in otu.index]
        tags = pd.read_csv(f"Raw_data/ok_108/tag_MF.csv", index_col=0)["Tag"]
        tags.index = [str(i) for i in tags.index]


    elif folder == "ok_30":
        otu = pd.read_csv(f"Data_for_git_process/ok_30/tax7_sum_relative.csv", index_col=0)
        otu.index = [str(i) for i in otu.index]
        tags = pd.read_csv(f"Data_for_git_process/ok_30/metadata_ok30_oshrit.csv", index_col=0)["high_lean"]
        tags.index = [str(i) for i in tags.index]
        folder_ = f"2D_otus_dendogram_ordered/{folder}"


    elif folder == "ok_94":
        otu = pd.read_csv(f"Data_for_git_process/ok_94/tax7_log_subpca.csv", index_col=0)
        otu.index = [str(i) for i in otu.index]
        tags = pd.read_csv(f"Data_for_git_process/ok_94/metadata_ok94_oshrit.tsv", index_col=0, sep="\t")["CD_or_UC"]
        tags.index = [str(i) for i in tags.index]
        tags[tags == "CD"] = 1.0
        tags[tags == "UC"] = 0.0
        tags[tags == "UC"] = 0.0
        tags = tags.astype(float)
        folder_ = f"2D_otus_dendogram_ordered/{folder}"


    elif folder == "GGMP":
        otu = pd.read_csv(f"Data_for_git_process/GGMP/tax7_relative_sum.csv", index_col=0)
        only_1000 = pd.read_csv("Raw_data/GGMP/only_1000.csv", index_col=0)
        common = list(only_1000.index.intersection(otu.index))
        otu = otu.loc[common]
        otu.index = [str(i) for i in otu.index]
        tags = pd.read_csv(f"Raw_data/GGMP/mapping_file.csv", index_col=0)["bmi"]
        tags.index = [i.split(".")[-1] for i in tags.index]
        tags = tags.loc[otu.index]



    elif folder == "PRJNA412501":
        otu = pd.read_csv(f"Data_for_git_process/PRJNA412501/tax7_relative_sum.csv", index_col=0)
        tags = pd.read_csv("Data_for_git_process/PRJNA412501/tag_inter.csv", index_col=0)["Tag"]
        tags.index = [f"{i}_1" for i in tags.index]
        folder_ = f"2D_otus_dendogram_ordered/{folder}"


    elif folder == "PRJNA419097":
        otu = pd.read_csv(f"Data_for_git_process/PRJNA419097/tax7_relative_sum.csv", index_col=0)
        tags = pd.read_csv("Data_for_git_process/PRJNA419097/tag_inter.csv", index_col=0)["Tag"]
        tags.index = [f"{i}_1" for i in tags.index]
        folder_ = f"2D_otus_dendogram_ordered/{folder}"


    elif folder == "PRJNA353587":
        otu = pd.read_csv(f"Data_for_git_process/PRJNA353587/tax7_relative_sum.csv", index_col=0)
        tags = pd.read_csv("Data_for_git_process/PRJNA353587/tag_inter.csv", index_col=0)["Tag"]
        folder_ = f"2D_otus_dendogram_ordered/{folder}"


    elif folder == "ERP020401":
        otu = pd.read_csv(f"Data_for_git_process/ERP020401/tax7_relative_sum.csv", index_col=0)
        tags = pd.read_csv("Data_for_git_process/ERP020401/final_meta.csv", index_col=0)["Tag"]
        folder_ = f"2D_otus_dendogram_ordered/{folder}"


    elif folder == "ERP021216":
        otu = pd.read_csv(f"Data_for_git_process/ERP021216/tax7_relative_sum.csv", index_col=0)
        tags = pd.read_csv("Data_for_git_process/ERP021216/final_meta.csv", index_col=0)["Tag"]
        folder_ = f"2D_otus_dendogram_ordered/{folder}"


    elif folder == "IBD/depth":
        otu = pd.read_csv(f"Data_for_git_process/{folder}/tax7_relative_sum_{p}.csv", index_col=0)
        otu.index = [str(i) for i in otu.index]
        tags = pd.read_csv("Data/IBD/tag_IBD_VS_ALL.csv", index_col=0)["Tag"]
        tags.index = [str(i) for i in tags.index]


    elif folder == "jacob/depth":
        otu = pd.read_csv(f"Data_for_git_process/{folder}/tax7_sum_relative_{p}.csv", index_col=0)
        tags = pd.read_csv(f"Data_for_git_process/jacob/tag.csv", index_col=0)


    elif folder == "ibd/remove_species":
        otu = pd.read_csv(f"Data_for_git_process/{folder}/tax7_sum_relative_{p}.csv", index_col=0)
        tags = pd.read_csv(f"Data/IBD/tag_IBD_VS_ALL.csv", index_col=0)["Tag"]
        # tags.index = [str(i) for i in tags.index]
        # imgs, names = load_img(f"2D_otus_dendogram_ordered/ibd_remove_species_{p}", tags)
        # tags = tags.reindex(names)

    elif folder == "jacob/remove_species":
        otu = pd.read_csv(f"Data_for_git_process/{folder}/tax7_relative_sum_{p}.csv", index_col=0)
        tags = pd.read_csv(f"Data_for_git_process/jacob/tag.csv", index_col=0)


    elif folder == "ok_94/remove_species":
        otu = pd.read_csv(f"Data_for_git_process/{folder}/tax7_relative_sum_{p}.csv", index_col=0)
        tags = pd.read_csv(f"Data_for_git_process/ok_94/metadata_ok94_oshrit.tsv", index_col=0, sep="\t")["CD_or_UC"]
        tags[tags == "CD"] = 1.0
        tags[tags == "UC"] = 0.0
        tags[tags == "UC"] = 0.0
        tags = tags.astype(float)

    elif folder == "he/remove_species":
        otu = pd.read_csv(f"Data_for_git_process/{folder}/tax7_relative_sum_{p}.csv", index_col=0)
        tags = pd.read_csv(f"Data_for_git_process/he/tag.csv", index_col=0)
        tags.index = ["12021." + i.replace("-", ".") for i in tags.index]


    elif folder == "Cirrhosis_Knight_Lab/remove_species":
        otu = pd.read_csv(f"Data_for_git_process/{folder}/tax7_relative_sum_{p}.csv", index_col=0)
        tags = pd.read_csv(f"Data/Cirrhosis_Knight_Lab/tag.csv", index_col=0)


    elif folder == "PRJNA353587/remove_species":
        otu = pd.read_csv(f"Data_for_git_process/{folder}/tax7_relative_sum_{p}.csv", index_col=0)
        tags = pd.read_csv("Data_for_git_process/PRJNA353587/tag_inter.csv", index_col=0)["Tag"]

    elif folder == "ERP020401/remove_species":
        otu = pd.read_csv(f"Data_for_git_process/{folder}/tax7_relative_sum_{p}.csv", index_col=0)
        tags = pd.read_csv("Data_for_git_process/ERP020401/final_meta.csv", index_col=0)["Tag"]


    elif folder == "ERP021216/remove_species":
        otu = pd.read_csv(f"Data_for_git_process/{folder}/tax7_relative_sum_{p}.csv", index_col=0)
        tags = pd.read_csv("Data_for_git_process/ERP021216/final_meta.csv", index_col=0)["Tag"]


    elif folder == "White_vs_black_vagina/remove_species":
        otu = pd.read_csv(f"Data_for_git_process/{folder}/tax7_relative_sum_{p}.csv", index_col=0)
        tags = pd.read_csv(f"Data/White_vs_black_vagina/tag.csv", index_col=0)


    elif folder == "Cirrhosis/remove_species":
        otu = pd.read_csv(f"Data_for_git_process/{folder}/tax7_relative_sum_{p}.csv", index_col=0)
        otu.index = [str(i) for i in otu.index]
        tags = pd.read_csv(f"Data/PopPhy_data/Cirrhosis/tag.csv", index_col=0)["Tag"]
        tags.index = [str(i) for i in tags.index]


    elif folder == "Male_vs_female/remove_species":
        otu = pd.read_csv(f"Data_for_git_process/{folder}/tax7_relative_sum_{p}.csv", index_col=0)
        tags = pd.read_csv(f"Data/Male_vs_female/tag.csv", index_col=0)


    elif folder == "PRJNA419097/remove_species":
        otu = pd.read_csv(f"Data_for_git_process/{folder}/tax7_relative_sum_{p}.csv", index_col=0)
        tags = pd.read_csv("Data_for_git_process/PRJNA419097/tag_inter.csv", index_col=0)["Tag"]
        tags.index = [f"{i}_1" for i in tags.index]


    elif folder == "data1_metab/remove_species":
        otu = pd.read_csv(f"Data_for_git_process/{folder}/tax7_relative_sum_{p}.csv", index_col=0)
        tags = pd.read_csv(f"Data_for_git_process/data1_metab/tag_for_learning.csv", index_col=0)["Tag"]


    elif folder == "ok_94/depth":
        otu = pd.read_csv(f"Data_for_git_process/{folder}/tax7_sum_relative_{p}.csv", index_col=0)
        tags = pd.read_csv(f"Data_for_git_process/ok_94/metadata_ok94_oshrit.tsv", index_col=0, sep="\t")["CD_or_UC"]
        tags[tags == "CD"] = 1.0
        tags[tags == "UC"] = 0.0
        tags[tags == "UC"] = 0.0
        tags = tags.astype(float)

    elif folder == "White_vs_black_vagina/depth":
        otu = pd.read_csv(f"Data_for_git_process/{folder}/tax7_sum_relative_{p}.csv", index_col=0)
        tags = pd.read_csv(f"Data/White_vs_black_vagina/tag.csv", index_col=0)


    elif folder == "Cirrhosis_Knight_Lab/depth":
        otu = pd.read_csv(f"Data_for_git_process/{folder}/tax7_sum_relative_{p}.csv", index_col=0)
        tags = pd.read_csv(f"Data/Cirrhosis_Knight_Lab/tag.csv", index_col=0)


    elif folder == "Male_vs_female/depth":
        otu = pd.read_csv(f"Data_for_git_process/{folder}/tax7_sum_relative_{p}.csv", index_col=0)
        tags = pd.read_csv(f"Data/Male_vs_female/tag.csv", index_col=0)


    elif folder == "Cirrhosis/depth":
        otu = pd.read_csv(f"Data_for_git_process/{folder}/tax7_sum_relative_{p}.csv", index_col=0)
        otu.index = [str(i) for i in otu.index]
        tags = pd.read_csv(f"Data/PopPhy_data/Cirrhosis/tag.csv", index_col=0)["Tag"]
        tags.index = [str(i) for i in tags.index]


    if shuffle:
        if folder == "GGMP" or folder == "Cirrhosis_Knight_Lab" or folder == "jacob" or folder == "he" or folder == "White_vs_black_vagina" or folder == "Male_vs_female" or folder == "data2_metab" or folder == "data3_metab" or folder == "data5_metab" or folder == "data6_metab" or folder == "T2D":
            for col in tags.columns:
                np.random.shuffle(tags[col].values)
        else:
            tags = tags.to_frame()
            for col in tags.columns:
                np.random.shuffle(tags[col].values)

    return tags, otu, folder_

def calc_sis_corr(significant,relative):
    my_list = [i.replace("s__","").replace("g__","").replace("f__","").replace("o__","").replace("c__","").replace("p__","").replace("k__","").split(";") for i in significant.index]
    cleaned_list_of_lists = [[item for item in inner_list if item != ""] for inner_list in my_list]
    for i in cleaned_list_of_lists:
        if len(i) == 1:
            cleaned_list_of_lists.remove(i)
    fathers = [i[-2] for i in cleaned_list_of_lists]


def fix_bact(df):
    df.index = [i.replace("k__","").replace("p__","").replace("c__","").replace("o__","").replace("f__","").replace("g__","").replace("s__","").replace("; ; ; ; ; ; ","").replace("; ; ; ; ; ","").replace("; ; ; ; ","").replace("; ; ; ","").replace("; ; ","").replace(";;;;;;","").replace(";;;;;","").replace(";;;; ","").replace(";;;;","").replace(";;;","").replace(";;","") for i in df.index]
    return df

def make_pairs(df):
    all_sccs = []
    for father, df_ in df.groupby("father"):
        if len(df_.index) >1:
            del df_["len"]
            del df_["father"]
            # Create pairs combinations from the values in column "0"
            pairs = list(itertools.combinations(df_.index, 2))
            # Create a new DataFrame with columns "s1" and "s2" from the pairs
            new_df = pd.DataFrame(index=df_.columns, columns=["s1", "s2"])
            for pair in pairs:
               new_df["s1"] = df_.loc[pair[0]]
               new_df["s2"] = df_.loc[pair[1]]
               scc = spearmanr(new_df["s1"],new_df["s2"])[0]
               all_sccs.append(scc)

    return all_sccs

def calc_father(df):
    if len(df.index)>0:
        df["father"] = [i.split(";")[-2] for i in df.index]
        return df
    else:
        return None


def find_sisters(significant,df,higher=None):
    df = df[significant.index]
    df =df.T
    df = fix_bact(df)
    my_list = [i.split(";") for i in df.index]
    cleaned_list_of_lists = [[item for item in inner_list if item != ""] for inner_list in my_list]
    df["len"] = [len(i) for i in cleaned_list_of_lists]
    df7 = df[df["len"] == 7]
    df6 = df[df["len"] == 6]
    df5 = df[df["len"] == 5]
    df4 = df[df["len"] == 4]
    df3 = df[df["len"] == 3]
    df2 = df[df["len"] == 2]

    try:
        df7 = calc_father(df7)
        sis7 = make_pairs(df7)
    except:
       sis7 = []
    try:
        df6 = calc_father(df6)
        sis6 = make_pairs(df6)
    except:
        sis6 = []
    try:
        df5= calc_father(df5)
        sis5 = make_pairs(df5)
    except:
        sis5 = []
    try:
        df4 =  calc_father(df4)
        sis4 = make_pairs(df4)
    except:
        sis4 = []
    try:
        df3 =  calc_father(df3)
        sis3 = make_pairs(df3)
    except:
        sis3 = []
    try:
        df2= calc_father(df2)
        sis2 = make_pairs(df2)
    except:
        sis2 = []

    all_corrs = list()
    if sis2:
        all_corrs.extend(sis2)
    if sis3:
        all_corrs.extend(sis3)
    if sis4:
        all_corrs.extend(sis4)
    if sis5:
        all_corrs.extend(sis5)
    if sis6:
        all_corrs.extend(sis6)
    if sis7:
        all_corrs.extend(sis7)
    return all_corrs

if __name__ == '__main__':
    mpl.rc('font', family='Times New Roman')
    plt.figure(figsize=(5, 5))
    SIZE = 25
    SHUFFLE = False
    all_corrs = list()
    for FOLDER in ["IBD","jacob","he","Cirrhosis","Cirrhosis_Knight_Lab","ok_30","ok_94","ERP020401","PRJNA419097","PRJNA353587","Obesity","T2D"]:
        print(FOLDER)
        df = pd.read_csv(f"corr_results/man/{FOLDER}.csv",index_col=0)
        significant = df[df["p"] < 0.05]
        meta, relative, folder = load_data(FOLDER, None, shuffle=SHUFFLE)
        corrs = find_sisters(significant,relative)
        all_corrs.append(corrs)
    flattend = list(chain(*all_corrs))
    plt.hist(flattend,bins=30,color="grey")
    MEAN = np.mean(flattend)
    plt.axvline(x=MEAN, color='deeppink', linestyle='dashed', linewidth=2, label=f'Mean = {MEAN:.2f}')
    plt.axvline(x=0, color='black', linewidth=2)
    plt.legend(fontsize=18)
    plt.xlabel("Sisters' SCC",fontsize=SIZE)
    plt.ylabel("Frequency",fontsize=SIZE)
    plt.xticks([-0.5,0,0.5,1], fontsize=SIZE)

    plt.yticks(fontsize=SIZE)
    plt.tight_layout()
    plt.savefig("corr_results/sis_distribution.png")
    plt.show()
    c=0
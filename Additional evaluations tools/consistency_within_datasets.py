import pandas as pd
def fix_ibd_taxa(df,model):
    if model == "deseq" or model == "ancom":
        df.index = [i.replace(";","; ") for i in df.index]

    elif model == "mimic":
        list_leaves = list(df.index)
        list_lens = [len(i.split(";")) for i in list_leaves]
        otu_train_cols = list()
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
        df.index = otu_train_cols
        df.index = [i.replace(";", "; ") for i in df.index]

    elif model == "aldex" or model == "ancom-bc2" or model == "ancom-bc2_shuffled":
        df.index = [i.replace(".", "..") for i in df.index]

    return df

def calc_percent_for_plot(df):
    num_s = len(df)
    counts = df['Total'].value_counts()
    counts = counts.to_frame('Counts')

    # Create a DataFrame with the desired index values
    index_values = [1, 2, 3, 4, 5]
    counts_df = pd.DataFrame(index=index_values)
    counts_df = pd.concat([counts_df, counts], axis=1)
    counts_df = counts_df.fillna(0.0)
    counts_df = (counts_df / num_s) * 100
    return counts_df

if __name__ == '__main__':
    DATA = "IBD"
    p = 0.05
    col = "we"#"wi"#"we"
    models_list = ["mimic_relative","mimic","linda_results","deseq","ancom","aldex","ancom-bc2","lefse"]
    MODEL = "ancom-bc2_shuffled"
    CORRECTED = False
    ibd_datasets_names = ["IBD", "ERP021216", "ok_94", "PRJNA353587", "PRJNA419097"]
    if MODEL == "linda_results":
        df2 =pd.read_csv(f"{MODEL}/LINDA_ERP021216.csv",index_col=0)
        df3 =pd.read_csv(f"{MODEL}/LINDA_ok_94.csv",index_col=0)
        df4 =pd.read_csv(f"{MODEL}/LINDA_PRJNA353587.csv",index_col=0)
        df5 =pd.read_csv(f"{MODEL}/LINDA_PRJNA419097.csv",index_col=0)
        common = list(df2.index.intersection(df3.index).intersection(df4.index).intersection(df5.index))
        df2 = df2.loc[common]
        df3 = df3.loc[common]
        df4 = df4.loc[common]
        df5 = df5.loc[common]
        if CORRECTED:
            p2 = df2["output.Tag.padj"]
            p3 = df3["output.CD_or_UC.padj"]
            p4 = df4["output.Tag.padj"]
            p5 = df5["output.Tag.padj"]
        else:
            p2 = df2["output.Tag.pvalue"]
            p3 = df3["output.CD_or_UC.pvalue"]
            p4 = df4["output.Tag.pvalue"]
            p5 = df5["output.Tag.pvalue"]

        # SIGNIFICANT
        s2 = p2[p2.values < p]
        s3 = p3[p3.values < p]
        s4 = p4[p4.values < p]
        s5 = p5[p5.values < p]
        significance_dist = pd.DataFrame(index=common,columns=['D2','D3','D4','D5'])
        for specie in common:
            if specie in s2.index:
                significance_dist['D2'][specie] = 1.0
            else:
                significance_dist['D2'][specie] = 0.0
            if specie in s3.index:
                significance_dist['D3'][specie] = 1.0
            else:
                significance_dist['D3'][specie] = 0.0
            if specie in s4.index:
                significance_dist['D4'][specie] = 1.0
            else:
                significance_dist['D4'][specie] = 0.0
            if specie in s5.index:
                significance_dist['D5'][specie] = 1.0
            else:
                significance_dist['D5'][specie] = 0.0

        significance_dist["Total"] = significance_dist.sum(axis=1)
        for_hist = significance_dist[significance_dist["Total"]!=0.0]
        if CORRECTED:
            for_hist.to_csv(f"rev_fig/consistency_of_model_between_datasets/{DATA}/linda-c.csv")
        else:
            for_hist.to_csv(f"rev_fig/consistency_of_model_between_datasets/{DATA}/linda.csv")

    elif MODEL == "deseq":
        if CORRECTED:
            p1 = pd.read_csv(f"ancom_bc2/data/IBD/{MODEL}_p_vals_corrected.csv",index_col=0)
            p2 = pd.read_csv(f"ancom_bc2/data/ERP021216/{MODEL}_p_vals_corrected.csv",index_col=0)
            p3 = pd.read_csv(f"ancom_bc2/data/ok_94/{MODEL}_p_vals_corrected.csv",index_col=0)
            p4 = pd.read_csv(f"ancom_bc2/data/PRJNA353587/{MODEL}_p_vals_corrected.csv",index_col=0)
            p5 = pd.read_csv(f"ancom_bc2/data/PRJNA419097/{MODEL}_p_vals_corrected.csv",index_col=0)
        else:
            p1 = pd.read_csv(f"ancom_bc2/data/IBD/{MODEL}_p_vals.csv", index_col=0)
            p2 = pd.read_csv(f"ancom_bc2/data/ERP021216/{MODEL}_p_vals.csv", index_col=0)
            p3 = pd.read_csv(f"ancom_bc2/data/ok_94/{MODEL}_p_vals.csv", index_col=0)
            p4 = pd.read_csv(f"ancom_bc2/data/PRJNA353587/{MODEL}_p_vals.csv", index_col=0)
            p5 = pd.read_csv(f"ancom_bc2/data/PRJNA419097/{MODEL}_p_vals.csv", index_col=0)

        p1 = fix_ibd_taxa(p1,MODEL)
        common = list(p1.index.intersection(p2.index).intersection(p3.index).intersection(p4.index).intersection(p5.index))
        p1 = p1.loc[common]
        p2 = p2.loc[common]
        p3 = p3.loc[common]
        p4 = p4.loc[common]
        p5 = p5.loc[common]
        # SIGNIFICANT
        s1 = p1[p1.values < p]
        s2 = p2[p2.values < p]
        s3 = p3[p3.values < p]
        s4 = p4[p4.values < p]
        s5 = p5[p5.values < p]
        significance_dist = pd.DataFrame(index=common, columns=['D1','D2', 'D3', 'D4', 'D5'])
        for specie in common:
            if specie in s1.index:
                significance_dist['D1'][specie] = 1.0
            else:
                significance_dist['D1'][specie] = 0.0
            if specie in s2.index:
                significance_dist['D2'][specie] = 1.0
            else:
                significance_dist['D2'][specie] = 0.0
            if specie in s3.index:
                significance_dist['D3'][specie] = 1.0
            else:
                significance_dist['D3'][specie] = 0.0
            if specie in s4.index:
                significance_dist['D4'][specie] = 1.0
            else:
                significance_dist['D4'][specie] = 0.0
            if specie in s5.index:
                significance_dist['D5'][specie] = 1.0
            else:
                significance_dist['D5'][specie] = 0.0

        significance_dist["Total"] = significance_dist.sum(axis=1)
        for_hist = significance_dist[significance_dist["Total"] != 0.0]
        if CORRECTED:
            for_hist.to_csv(f"rev_fig/consistency_of_model_between_datasets/{DATA}/deseq-c.csv")
        else:
            for_hist.to_csv(f"rev_fig/consistency_of_model_between_datasets/{DATA}/deseq.csv")

    elif MODEL == "ancom":
        if CORRECTED:
            p1 = pd.read_csv(f"ancom_bc2/data/IBD/{MODEL}_p_vals_corrected.csv",index_col=0)
            p2 = pd.read_csv(f"ancom_bc2/data/ERP021216/{MODEL}_p_vals_corrected.csv",index_col=0)
            p3 = pd.read_csv(f"ancom_bc2/data/ok_94/{MODEL}_p_vals_corrected.csv",index_col=0)
            p4 = pd.read_csv(f"ancom_bc2/data/PRJNA353587/{MODEL}_p_vals_corrected.csv",index_col=0)
            p5 = pd.read_csv(f"ancom_bc2/data/PRJNA419097/{MODEL}_p_vals_corrected.csv",index_col=0)
        else:
            p1 = pd.read_csv(f"ancom_bc2/data/IBD/{MODEL}_p_vals.csv", index_col=0)
            p2 = pd.read_csv(f"ancom_bc2/data/ERP021216/{MODEL}_p_vals.csv", index_col=0)
            p3 = pd.read_csv(f"ancom_bc2/data/ok_94/{MODEL}_p_vals.csv", index_col=0)
            p4 = pd.read_csv(f"ancom_bc2/data/PRJNA353587/{MODEL}_p_vals.csv", index_col=0)
            p5 = pd.read_csv(f"ancom_bc2/data/PRJNA419097/{MODEL}_p_vals.csv", index_col=0)

        p1 = fix_ibd_taxa(p1,MODEL)
        common = list(
            p1.index.intersection(p2.index).intersection(p3.index).intersection(p4.index).intersection(p5.index))
        # SIGNIFICANT
        s1 = p1[p1.values == True]
        s2 = p2[p2.values == True]
        s3 = p3[p3.values == True]
        s4 = p4[p4.values == True]
        s5 = p5[p5.values == True]
        significance_dist = pd.DataFrame(index=common, columns=['D1', 'D2', 'D3', 'D4', 'D5'])
        for specie in common:
            if specie in s1.index:
                significance_dist['D1'][specie] = 1.0
            else:
                significance_dist['D1'][specie] = 0.0
            if specie in s2.index:
                significance_dist['D2'][specie] = 1.0
            else:
                significance_dist['D2'][specie] = 0.0
            if specie in s3.index:
                significance_dist['D3'][specie] = 1.0
            else:
                significance_dist['D3'][specie] = 0.0
            if specie in s4.index:
                significance_dist['D4'][specie] = 1.0
            else:
                significance_dist['D4'][specie] = 0.0
            if specie in s5.index:
                significance_dist['D5'][specie] = 1.0
            else:
                significance_dist['D5'][specie] = 0.0

        significance_dist["Total"] = significance_dist.sum(axis=1)
        for_hist = significance_dist[significance_dist["Total"] != 0.0]
        if CORRECTED:
            for_hist.to_csv(f"rev_fig/consistency_of_model_between_datasets/{DATA}/{MODEL}-c.csv")
        else:
            for_hist.to_csv(f"rev_fig/consistency_of_model_between_datasets/{DATA}/{MODEL}.csv")

    elif MODEL == "aldex":
        if CORRECTED:
            p1 = pd.read_csv(f"ancom_bc2/data/IBD/{MODEL}_real_holm.csv",index_col=0)[f"{col}.eBH"]
            p2 = pd.read_csv(f"ancom_bc2/data/ERP021216/{MODEL}_real_holm.csv",index_col=0)[f"{col}.eBH"]
            p3 = pd.read_csv(f"ancom_bc2/data/ok_94/{MODEL}_real_holm.csv",index_col=0)[f"{col}.eBH"]
            p4 = pd.read_csv(f"ancom_bc2/data/PRJNA353587/{MODEL}_real_holm.csv",index_col=0)[f"{col}.eBH"]
            p5 = pd.read_csv(f"ancom_bc2/data/PRJNA419097/{MODEL}_real_holm.csv",index_col=0)[f"{col}.eBH"]
        else:
            p1 = pd.read_csv(f"ancom_bc2/data/IBD/{MODEL}_real_holm.csv", index_col=0)[f"{col}.ep"]
            p2 = pd.read_csv(f"ancom_bc2/data/ERP021216/{MODEL}_real_holm.csv", index_col=0)[f"{col}.ep"]
            p3 = pd.read_csv(f"ancom_bc2/data/ok_94/{MODEL}_real_holm.csv", index_col=0)[f"{col}.ep"]
            p4 = pd.read_csv(f"ancom_bc2/data/PRJNA353587/{MODEL}_real_holm.csv", index_col=0)[f"{col}.ep"]
            p5 = pd.read_csv(f"ancom_bc2/data/PRJNA419097/{MODEL}_real_holm.csv", index_col=0)[f"{col}.ep"]

        p1 = fix_ibd_taxa(p1,MODEL)
        common = list(
            p1.index.intersection(p2.index).intersection(p3.index).intersection(p4.index).intersection(p5.index))
        # SIGNIFICANT
        s1 = p1[p1.values < p]
        s2 = p2[p2.values < p]
        s3 = p3[p3.values < p]
        s4 = p4[p4.values < p]
        s5 = p5[p5.values < p]
        significance_dist = pd.DataFrame(index=common, columns=['D1', 'D2', 'D3', 'D4', 'D5'])
        for specie in common:
            if specie in s1.index:
                significance_dist['D1'][specie] = 1.0
            else:
                significance_dist['D1'][specie] = 0.0
            if specie in s2.index:
                significance_dist['D2'][specie] = 1.0
            else:
                significance_dist['D2'][specie] = 0.0
            if specie in s3.index:
                significance_dist['D3'][specie] = 1.0
            else:
                significance_dist['D3'][specie] = 0.0
            if specie in s4.index:
                significance_dist['D4'][specie] = 1.0
            else:
                significance_dist['D4'][specie] = 0.0
            if specie in s5.index:
                significance_dist['D5'][specie] = 1.0
            else:
                significance_dist['D5'][specie] = 0.0

        significance_dist["Total"] = significance_dist.sum(axis=1)
        for_hist = significance_dist[significance_dist["Total"] != 0.0]
        if CORRECTED:
            for_hist.to_csv(f"rev_fig/consistency_of_model_between_datasets/{DATA}/{MODEL}_{col}-c.csv")
        else:
            for_hist.to_csv(f"rev_fig/consistency_of_model_between_datasets/{DATA}/{MODEL}_{col}.csv")

    elif MODEL == "ancom-bc2":
        if CORRECTED:
            p1 = pd.read_csv(f"ancom_bc2/data/IBD/real_holm.csv", index_col=0)["q_Tag"]
            p2 = pd.read_csv(f"ancom_bc2/data/ERP021216/real_holm.csv", index_col=0)["q_Tag"]
            p3 = pd.read_csv(f"ancom_bc2/data/ok_94/real_holm.csv", index_col=0)["q_Tag"]
            p4 = pd.read_csv(f"ancom_bc2/data/PRJNA353587/real_holm.csv", index_col=0)["q_Tag"]
            p5 = pd.read_csv(f"ancom_bc2/data/PRJNA419097/real_holm.csv", index_col=0)["q_Tag"]
        else:
            p1 = pd.read_csv(f"ancom_bc2/data/IBD/real_holm.csv", index_col=0)["p_Tag"]
            p2 = pd.read_csv(f"ancom_bc2/data/ERP021216/real_holm.csv", index_col=0)["p_Tag"]
            p3 = pd.read_csv(f"ancom_bc2/data/ok_94/real_holm.csv", index_col=0)["p_Tag"]
            p4 = pd.read_csv(f"ancom_bc2/data/PRJNA353587/real_holm.csv", index_col=0)["p_Tag"]
            p5 = pd.read_csv(f"ancom_bc2/data/PRJNA419097/real_holm.csv", index_col=0)["p_Tag"]

        p1 = fix_ibd_taxa(p1, MODEL)
        common = list(
            p1.index.intersection(p2.index).intersection(p3.index).intersection(p4.index).intersection(p5.index))
        # SIGNIFICANT
        s1 = p1[p1.values < p]
        s2 = p2[p2.values < p]
        s3 = p3[p3.values < p]
        s4 = p4[p4.values < p]
        s5 = p5[p5.values < p]
        significance_dist = pd.DataFrame(index=common, columns=['D1', 'D2', 'D3', 'D4', 'D5'])
        for specie in common:
            if specie in s1.index:
                significance_dist['D1'][specie] = 1.0
            else:
                significance_dist['D1'][specie] = 0.0
            if specie in s2.index:
                significance_dist['D2'][specie] = 1.0
            else:
                significance_dist['D2'][specie] = 0.0
            if specie in s3.index:
                significance_dist['D3'][specie] = 1.0
            else:
                significance_dist['D3'][specie] = 0.0
            if specie in s4.index:
                significance_dist['D4'][specie] = 1.0
            else:
                significance_dist['D4'][specie] = 0.0
            if specie in s5.index:
                significance_dist['D5'][specie] = 1.0
            else:
                significance_dist['D5'][specie] = 0.0

        significance_dist["Total"] = significance_dist.sum(axis=1)
        for_hist = significance_dist[significance_dist["Total"] != 0.0]
        if CORRECTED:
            for_hist.to_csv(f"rev_fig/consistency_of_model_between_datasets/{DATA}/{MODEL}-c.csv")
        else:
            for_hist.to_csv(f"rev_fig/consistency_of_model_between_datasets/{DATA}/{MODEL}.csv")

    elif MODEL == "ancom-bc2_shuffled":
        if CORRECTED:
            p1 = pd.read_csv(f"ancom_bc2/data/IBD/shuffle_holm.csv", index_col=0)["q_Tag"]
            p2 = pd.read_csv(f"ancom_bc2/data/ERP021216/shuffle_holm.csv", index_col=0)["q_Tag"]
            p3 = pd.read_csv(f"ancom_bc2/data/ok_94/shuffle_holm.csv", index_col=0)["q_Tag"]
            p4 = pd.read_csv(f"ancom_bc2/data/PRJNA353587/shuffle_holm.csv", index_col=0)["q_Tag"]
            p5 = pd.read_csv(f"ancom_bc2/data/PRJNA419097/shuffle_holm.csv", index_col=0)["q_Tag"]
        else:
            p1 = pd.read_csv(f"ancom_bc2/data/IBD/shuffle_holm.csv", index_col=0)["p_Tag"]
            p2 = pd.read_csv(f"ancom_bc2/data/ERP021216/shuffle_holm.csv", index_col=0)["p_Tag"]
            p3 = pd.read_csv(f"ancom_bc2/data/ok_94/shuffle_holm.csv", index_col=0)["p_Tag"]
            p4 = pd.read_csv(f"ancom_bc2/data/PRJNA353587/shuffle_holm.csv", index_col=0)["p_Tag"]
            p5 = pd.read_csv(f"ancom_bc2/data/PRJNA419097/shuffle_holm.csv", index_col=0)["p_Tag"]

        p1 = fix_ibd_taxa(p1, MODEL)
        common = list(
            p1.index.intersection(p2.index).intersection(p3.index).intersection(p4.index).intersection(p5.index))
        # SIGNIFICANT
        s1 = p1[p1.values < p]
        s2 = p2[p2.values < p]
        s3 = p3[p3.values < p]
        s4 = p4[p4.values < p]
        s5 = p5[p5.values < p]

        significance_dist = pd.DataFrame(index=common, columns=['D1', 'D2', 'D3', 'D4', 'D5'])
        for specie in common:
            if specie in s1.index:
                significance_dist['D1'][specie] = 1.0
            else:
                significance_dist['D1'][specie] = 0.0
            if specie in s2.index:
                significance_dist['D2'][specie] = 1.0
            else:
                significance_dist['D2'][specie] = 0.0
            if specie in s3.index:
                significance_dist['D3'][specie] = 1.0
            else:
                significance_dist['D3'][specie] = 0.0
            if specie in s4.index:
                significance_dist['D4'][specie] = 1.0
            else:
                significance_dist['D4'][specie] = 0.0
            if specie in s5.index:
                significance_dist['D5'][specie] = 1.0
            else:
                significance_dist['D5'][specie] = 0.0

        significance_dist["Total"] = significance_dist.sum(axis=1)
        for_hist = significance_dist[significance_dist["Total"] != 0.0]
        if CORRECTED:
            for_hist.to_csv(f"rev_fig/consistency_of_model_between_datasets/{DATA}/{MODEL}_shuffled-c.csv")
        else:
            for_hist.to_csv(f"rev_fig/consistency_of_model_between_datasets/{DATA}/{MODEL}_shuffled.csv")

    elif MODEL == "lefse":
        p1 = pd.read_csv(f"ancom_bc2/data/IBD/IBD_TRUE_lefse.csv", index_col=0)["P"]
        p2 = pd.read_csv(f"ancom_bc2/data/ERP021216/ERP021216_True_lefse.csv", index_col=0)["P"]
        p3 = pd.read_csv(f"ancom_bc2/data/ok_94/ok_94_True_lefse.csv", index_col=0)["P"]
        p4 = pd.read_csv(f"ancom_bc2/data/PRJNA353587/PRJNA353587_True_lefse.csv", index_col=0)["P"]
        p5 = pd.read_csv(f"ancom_bc2/data/PRJNA419097/PRJNA419097_True_lefse.csv", index_col=0)["P"]

        common = list(
            p1.index.intersection(p2.index).intersection(p3.index).intersection(p4.index).intersection(p5.index))
        # SIGNIFICANT
        s1 = p1[p1.values != "-"]
        s2 = p2[p2.values != "-"]
        s3 = p3[p3.values != "-"]
        s4 = p4[p4.values != "-"]
        s5 = p5[p5.values != "-"]
        significance_dist = pd.DataFrame(index=common, columns=['D1', 'D2', 'D3', 'D4', 'D5'])
        for specie in common:
            if specie in s1.index:
                significance_dist['D1'][specie] = 1.0
            else:
                significance_dist['D1'][specie] = 0.0
            if specie in s2.index:
                significance_dist['D2'][specie] = 1.0
            else:
                significance_dist['D2'][specie] = 0.0
            if specie in s3.index:
                significance_dist['D3'][specie] = 1.0
            else:
                significance_dist['D3'][specie] = 0.0
            if specie in s4.index:
                significance_dist['D4'][specie] = 1.0
            else:
                significance_dist['D4'][specie] = 0.0
            if specie in s5.index:
                significance_dist['D5'][specie] = 1.0
            else:
                significance_dist['D5'][specie] = 0.0

        significance_dist["Total"] = significance_dist.sum(axis=1)
        for_hist = significance_dist[significance_dist["Total"] != 0.0]

    elif MODEL == "mimic":

        # SIGNIFICANT
        s1 = pd.read_csv(f"ancom_bc2/data/IBD/df_corrs.csv", index_col=0)["p"]
        s2 = pd.read_csv(f"ancom_bc2/data/ERP021216/mimic_df_corrs.csv", index_col=0)["p"]
        s3 = pd.read_csv(f"ancom_bc2/data/ok_94/mimic_df_corrs.csv", index_col=0)["p"]
        s4 = pd.read_csv(f"ancom_bc2/data/PRJNA353587/mimic_df_corrs.csv", index_col=0)["p"]
        s5 = pd.read_csv(f"ancom_bc2/data/PRJNA419097/mimic_df_corrs.csv", index_col=0)["p"]

        s1 = fix_ibd_taxa(s1,"mimic")
        s2 = fix_ibd_taxa(s2,"mimic")
        s3 = fix_ibd_taxa(s3,"mimic")
        s4 = fix_ibd_taxa(s4,"mimic")
        s5 = fix_ibd_taxa(s5,"mimic")

        union_list = set(s1.index) | set(s2.index) | set(s3.index) | set(s4.index) | set(s5.index)
        common = list(union_list)

        significance_dist = pd.DataFrame(index=common, columns=['D1', 'D2', 'D3', 'D4', 'D5'])
        for specie in common:
            if specie in s1.index:
                significance_dist['D1'][specie] = 1.0
            else:
                significance_dist['D1'][specie] = 0.0
            if specie in s2.index:
                significance_dist['D2'][specie] = 1.0
            else:
                significance_dist['D2'][specie] = 0.0
            if specie in s3.index:
                significance_dist['D3'][specie] = 1.0
            else:
                significance_dist['D3'][specie] = 0.0
            if specie in s4.index:
                significance_dist['D4'][specie] = 1.0
            else:
                significance_dist['D4'][specie] = 0.0
            if specie in s5.index:
                significance_dist['D5'][specie] = 1.0
            else:
                significance_dist['D5'][specie] = 0.0

        significance_dist["Total"] = significance_dist.sum(axis=1)
        for_hist = significance_dist[significance_dist["Total"] != 0.0]

    elif MODEL == "mimic_relative":

        # SIGNIFICANT
        s1 = pd.read_csv(f"ancom_bc2/data/IBD/mimic_df_corrs_relative.csv", index_col=0)["p"]
        s2 = pd.read_csv(f"ancom_bc2/data/ERP021216/mimic_df_corrs_relative.csv", index_col=0)["p"]
        s3 = pd.read_csv(f"ancom_bc2/data/ok_94/mimic_df_corrs_relative.csv", index_col=0)["p"]


        s1 = fix_ibd_taxa(s1,"mimic")
        s2 = fix_ibd_taxa(s2,"mimic")
        s3 = fix_ibd_taxa(s3,"mimic")


        union_list = set(s1.index) | set(s2.index) | set(s3.index)
        common = list(union_list)

        significance_dist = pd.DataFrame(index=common, columns=['D1', 'D2', 'D3', 'D4', 'D5'])
        for specie in common:
            if specie in s1.index:
                significance_dist['D1'][specie] = 1.0
            else:
                significance_dist['D1'][specie] = 0.0
            if specie in s2.index:
                significance_dist['D2'][specie] = 1.0
            else:
                significance_dist['D2'][specie] = 0.0
            if specie in s3.index:
                significance_dist['D3'][specie] = 1.0
            else:
                significance_dist['D3'][specie] = 0.0
            # if specie in s4.index:
            #     significance_dist['D4'][specie] = 1.0
            # else:
            #     significance_dist['D4'][specie] = 0.0
            # if specie in s5.index:
            #     significance_dist['D5'][specie] = 1.0
            # else:
            #     significance_dist['D5'][specie] = 0.0

        significance_dist["Total"] = significance_dist.sum(axis=1)
        for_hist = significance_dist[significance_dist["Total"] != 0.0]

    for_plot = calc_percent_for_plot(for_hist)
    if CORRECTED:
        for_plot.to_csv(f"rev_fig/consistency_of_model_between_datasets/IBD/real/{MODEL}-{col}-c.csv")
    else:
        for_plot.to_csv(f"rev_fig/consistency_of_model_between_datasets/IBD/real/{MODEL}-{col}.csv")




        print()


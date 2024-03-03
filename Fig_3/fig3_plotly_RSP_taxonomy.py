import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np


if __name__ == '__main__':
    SIZE = 15
    WGS = True
    mpl.rc('font', family='Times New Roman')
    fig, ax = plt.subplots(1, 1, figsize=[4, 4])
    list_16s = ["Jacob", "He", "IBD", "Cirrhosis", "MF", "ERP020401", "ERP021216", "OK30", "OK94", "PRJNA353587",
                "PRJNA419097", "WB"]
    list_wgs = ["Cirrhosis-pop", "T2D", "Obesity", "WGS1", "WGS2", "WGS3", "WGS4", "WGS5"]
    list_dfs = list()
    if WGS:
        LIST = list_wgs
    else:
        LIST = list_16s
    for NAME in LIST:
        df = pd.read_csv(f"temp/taxonomy/{NAME}.csv", index_col=0)
        list_dfs.append(df)

    stacked_df = pd.concat(list_dfs)

    # Calculate the mean along axis=0 (rows) while ignoring NaN values
    averaged_df = stacked_df.groupby(level=0).mean()
    averaged_df = averaged_df.T
    if WGS:
        std_df = (stacked_df.groupby(level=0).std()) / (np.sqrt(8))
    else:
        std_df = (stacked_df.groupby(level=0).std()) / (np.sqrt(12))

    std_df = std_df.T
    ALPHA = 0.1

    ax.plot(averaged_df.index, averaged_df[1], label='1', color='deeppink',linewidth=3.0)
    ax.fill_between(list(averaged_df.index), (averaged_df[1] - std_df[1]).values.astype(float),
                    (averaged_df[1] + std_df[1]).values.astype(float), alpha=ALPHA, color='deeppink')
    ax.plot(averaged_df.index, averaged_df[2], label='2', color='deeppink',linestyle='--',linewidth=3.0)
    ax.fill_between(list(averaged_df.index), (averaged_df[2] - std_df[2]).values.astype(float),
                    (averaged_df[2] + std_df[2]).values.astype(float), alpha=ALPHA, color='deeppink')
    ax.plot(averaged_df.index, averaged_df[3], label='3', color='deeppink',linestyle=':',linewidth=3.0)
    ax.fill_between(list(averaged_df.index), (averaged_df[3] - std_df[3]).values.astype(float),
                    (averaged_df[3] + std_df[3]).values.astype(float), alpha=ALPHA, color='deeppink')
    ax.plot(averaged_df.index, averaged_df[4], label='4', color='deeppink',linestyle='-.')
    ax.fill_between(list(averaged_df.index), (averaged_df[4] - std_df[4]).values.astype(float),
                    (averaged_df[4] + std_df[4]).values.astype(float), alpha=ALPHA, color='deeppink')
    ax.plot(averaged_df.index, averaged_df[5], label='5', color='deeppink',linestyle='--')
    ax.fill_between(list(averaged_df.index), (averaged_df[5] - std_df[5]).values.astype(float),
                    (averaged_df[5] + std_df[5]).values.astype(float), alpha=ALPHA, color='deeppink')
    ax.plot(averaged_df.index, averaged_df[6], label='6', color='deeppink',linestyle=':')
    ax.fill_between(list(averaged_df.index), (averaged_df[6] - std_df[6]).values.astype(float),
                    (averaged_df[6] + std_df[6]).values.astype(float), alpha=ALPHA, color='deeppink')
    ax.plot(averaged_df.index, averaged_df[7], label='7', color='hotpink')
    ax.fill_between(list(averaged_df.index), (averaged_df[7] - std_df[7]).values.astype(float),
                    (averaged_df[7] + std_df[7]).values.astype(float), alpha=ALPHA, color='hotpink')

    plt.legend()
    plt.xlabel(r'$\beta$', fontsize=SIZE)
    plt.ylabel(r"RSP($\beta$) score", fontsize=SIZE)
    x_ticks = averaged_df.index[::3]  # Show ticks in jumps of 3
    ax.set_xticks(x_ticks)
    ax.set_xticklabels(x_ticks, fontsize=SIZE)
    plt.yticks(fontsize=SIZE)

    if WGS:
        axins = ax.inset_axes([9.2, -0.58, 6, 0.60], transform=ax.transData)
    else:
        axins = ax.inset_axes([9.2, -0.54, 6, 0.60], transform=ax.transData)

    if WGS:
        hist = pd.read_csv("temp/outline_figures/Fig_3/number_datasets_wgs.csv", index_col=0)
    else:
        hist = pd.read_csv("temp/outline_figures/Fig_3/number_datasets_16.csv", index_col=0)
    axins.bar(hist.index, hist["Number of datasets"], color=["hotpink"])
    axins.set_xlabel("Taxonomy")
    axins.set_ylabel("Datasets")
    axins.set_xticks([1, 2, 3, 4, 5, 6, 7])
    if WGS:
        plt.title("WGS", fontsize=SIZE, fontweight="bold")
    else:
        plt.title("16S", fontsize=SIZE, fontweight="bold")

    plt.tight_layout()
    if WGS:
        plt.savefig("temp/outline_figures/Fig_3/wgs_taxonomy.png")
    else:
        plt.savefig("temp/outline_figures/Fig_3/16s_taxonomy.png")

    plt.show()

    c = 0

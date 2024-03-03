import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np

if __name__ == '__main__':
    WGS = True
    LEGEND = False
    SIZE=15
    plt.figure(figsize=(4,4))
    mpl.rc('font', family='Times New Roman')
    list_alphas = [0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9,
                   0.95, 1.0]

    list_models = ["DeSeq", "DeSeq-C", "ANCOM", "ANCOM-C", "LINDA", "LINDA-C","LEfSe","ANCOM-BC 2","ANCOM-BC 2-C","ALDEx (Welch)","ALDEx (Welch)-C","ALDEx (Wilcoxon)","ALDEx (Wilcoxon)-C"
                   ,"miMic","miMic relative","ada-ANCOM"
                   ]#"miMic-fs-first"
    for_heatmap = pd.DataFrame(index=list_alphas, columns=list_models)
    for_std = pd.DataFrame(index=list_alphas, columns=list_models)
    for alpha in list_alphas:

        df = pd.read_csv(f"temp/Alpha/rev_1m/rev_comparisons_{alpha}.csv"
                         , index_col=0)
        if WGS:
            df = df.loc[["WGS-1", "WGS-2", "WGS-3", "WGS-4", "WGS-5", "T2D", "Obesity", "Cirrhosis-POP"]]
        else:
            df = df.loc[
                ["Cirrhosis", "ERP020401", "ERP021216", "He", "IBD", "Jacob", "MF", "OK30", "OK94", "PRJNA353587",
                 "PRJNA419097", "WB"]]

        to_plot = df[list_models].mean()
        for_heatmap.loc[alpha] = to_plot
        to_std = df[list_models].std() / (np.sqrt(len(df.index)))
        for_std.loc[alpha] = to_std

    plt.plot(for_heatmap.index, for_heatmap["DeSeq"], label='DeSeq', color='salmon')
    plt.fill_between(list(for_heatmap.index), (for_heatmap["DeSeq"] - for_std["DeSeq"]).values.astype(float),
                     (for_heatmap["DeSeq"] + for_std["DeSeq"]).values.astype(float), alpha=0.3, color='salmon')
    plt.plot(for_heatmap.index, for_heatmap["DeSeq-C"], label='DeSeq-C', color='tomato',linestyle='--')
    plt.fill_between(list(for_heatmap.index), (for_heatmap["DeSeq-C"] - for_std["DeSeq-C"]).values.astype(float),
                     (for_heatmap["DeSeq-C"] + for_std["DeSeq-C"]).values.astype(float), alpha=0.3, color='tomato',
                     )
    plt.plot(for_heatmap.index, for_heatmap["LEfSe"], label='LEfSe', color='gold')
    plt.fill_between(list(for_heatmap.index), (for_heatmap["LEfSe"] - for_std["LEfSe"]).values.astype(float),
                     (for_heatmap["LEfSe"] + for_std["LEfSe"]).values.astype(float), alpha=0.3, color='gold',
                     )

    plt.plot(for_heatmap.index, for_heatmap["ANCOM"], label='ANCOM', color='mediumaquamarine')
    plt.fill_between(list(for_heatmap.index), (for_heatmap["ANCOM"] - for_std["ANCOM"]).values.astype(float),
                     (for_heatmap["ANCOM"] + for_std["ANCOM"]).values.astype(float), alpha=0.3, color='mediumaquamarine',
                     )
    plt.plot(for_heatmap.index, for_heatmap["ANCOM-C"], label='ANCOM-C', color='seagreen',linestyle='--')
    plt.fill_between(list(for_heatmap.index), (for_heatmap["ANCOM-C"] - for_std["ANCOM-C"]).values.astype(float),
                     (for_heatmap["ANCOM-C"] + for_std["ANCOM-C"]).values.astype(float), alpha=0.3,
                     color='seagreen',
                     )

    plt.plot(for_heatmap.index, for_heatmap["ANCOM-BC 2"], label='ANCOM-BC2', color='green')
    plt.fill_between(list(for_heatmap.index), (for_heatmap["ANCOM-BC 2"] - for_std["ANCOM-BC 2"]).values.astype(float),
                     (for_heatmap["ANCOM-BC 2"] + for_std["ANCOM-BC 2"]).values.astype(float), alpha=0.3,
                     color='green',
                     )

    plt.plot(for_heatmap.index, for_heatmap["ANCOM-BC 2-C"], label='ANCOM-BC2-C', color='darkgreen',linestyle='--')
    plt.fill_between(list(for_heatmap.index), (for_heatmap["ANCOM-BC 2-C"] - for_std["ANCOM-BC 2-C"]).values.astype(float),
                     (for_heatmap["ANCOM-BC 2-C"] + for_std["ANCOM-BC 2-C"]).values.astype(float), alpha=0.3,
                     color='darkgreen',
                     )

    plt.plot(for_heatmap.index, for_heatmap["LINDA"], label='LINDA', color='lightskyblue')
    plt.fill_between(list(for_heatmap.index), (for_heatmap["LINDA"] - for_std["LINDA"]).values.astype(float),
                     (for_heatmap["LINDA"] + for_std["LINDA"]).values.astype(float), alpha=0.3,color='lightskyblue',
                     )
    plt.plot(for_heatmap.index, for_heatmap["LINDA-C"], label='LINDA-C', color='cornflowerblue',linestyle='--')
    plt.fill_between(list(for_heatmap.index), (for_heatmap["LINDA-C"] - for_std["LINDA-C"]).values.astype(float),
                     (for_heatmap["LINDA-C"] + for_std["LINDA-C"]).values.astype(float), alpha=0.3,
                     color='cornflowerblue',
                     )

    plt.plot(for_heatmap.index, for_heatmap["ALDEx (Welch)"], label='ALDEx (Welch)', color='grey')
    plt.fill_between(list(for_heatmap.index), (for_heatmap["ALDEx (Welch)"] - for_std["ALDEx (Welch)"]).values.astype(float),
                     (for_heatmap["ALDEx (Welch)"] + for_std["ALDEx (Welch)"]).values.astype(float), alpha=0.3,
                     color='grey',
                     )

    plt.plot(for_heatmap.index, for_heatmap["ALDEx (Welch)-C"], label='ALDEx (Welch)-C', color='dimgrey',linestyle='--')
    plt.fill_between(list(for_heatmap.index),
                     (for_heatmap["ALDEx (Welch)-C"] - for_std["ALDEx (Welch)-C"]).values.astype(float),
                     (for_heatmap["ALDEx (Welch)-C"] + for_std["ALDEx (Welch)-C"]).values.astype(float), alpha=0.3,
                     color='dimgrey',
                     )

    plt.plot(for_heatmap.index, for_heatmap["ALDEx (Wilcoxon)"], label='ALDEx (Wilcoxon)', color='lightslategray',linestyle=':')
    plt.fill_between(list(for_heatmap.index),
                     (for_heatmap["ALDEx (Wilcoxon)"] - for_std["ALDEx (Wilcoxon)"]).values.astype(float),
                     (for_heatmap["ALDEx (Wilcoxon)"] + for_std["ALDEx (Wilcoxon)"]).values.astype(float), alpha=0.3,
                     color='lightslategray',
                     )

    plt.plot(for_heatmap.index, for_heatmap["ALDEx (Wilcoxon)-C"], label='ALDEx (Wilcoxon)-C', color='slategray',linestyle='-')
    plt.fill_between(list(for_heatmap.index),
                     (for_heatmap["ALDEx (Wilcoxon)-C"] - for_std["ALDEx (Wilcoxon)-C"]).values.astype(float),
                     (for_heatmap["ALDEx (Wilcoxon)-C"] + for_std["ALDEx (Wilcoxon)-C"]).values.astype(float), alpha=0.3,
                     color='slategray',
                     )

    plt.plot(for_heatmap.index, for_heatmap["miMic"], label='miMic', color='hotpink')
    plt.fill_between(list(for_heatmap.index), (for_heatmap["miMic"] - for_std["miMic"]).values.astype(float),
                     (for_heatmap["miMic"] + for_std["miMic"]).values.astype(float), alpha=0.3,
                     color='hotpink',
                     )#miMic-fs-first

    plt.plot(for_heatmap.index, for_heatmap["miMic relative"], label='miMic relative', color='mediumorchid')
    plt.fill_between(list(for_heatmap.index),
                     (for_heatmap["miMic relative"] - for_std["miMic relative"]).values.astype(float),
                     (for_heatmap["miMic relative"] + for_std["miMic relative"]).values.astype(float), alpha=0.3,
                     color='mediumorchid',
                     )
    if not WGS:
        plt.plot(for_heatmap.index, for_heatmap["ada-ANCOM"], label='ada-ANCOM', color='saddlebrown')
        plt.fill_between(list(for_heatmap.index),
                         (for_heatmap["ada-ANCOM"] - for_std["ada-ANCOM"]).values.astype(float),
                         (for_heatmap["ada-ANCOM"] + for_std["ada-ANCOM"]).values.astype(float), alpha=0.3,
                         color='saddlebrown',
                         )
    if LEGEND:
        # Hide the plot


        # Add legend and set its properties
        plt.legend(loc='center', fontsize=15, frameon=False)
        plt.axis('off')

        # Adjust the layout
        plt.tight_layout()

        # Show the legend
        plt.show()
    else:


        plt.xlabel(r'$\beta$',fontsize=SIZE)
        plt.ylabel(r"RSP($\beta$) score",fontsize=SIZE)
        plt.xticks(fontsize=SIZE)
        plt.yticks(fontsize=SIZE)
        if WGS:
            plt.title("WGS",fontsize=SIZE,fontweight="bold")
        else:
            plt.title("16S",fontsize=SIZE,fontweight="bold")

        plt.tight_layout()
        if WGS:
            plt.savefig("temp/outline_figures/rev1/Fig_3/wgs_comp_.png")
        else:
            plt.savefig("temp/outline_figures/rev1/Fig_3/16s_comp_.png")

        plt.show()



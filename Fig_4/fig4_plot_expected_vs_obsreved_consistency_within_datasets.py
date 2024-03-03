import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
def plot_real_exp(real,exp,MODEL,color,avg=False):
    SIZE = 20
    if not avg:
        both = pd.DataFrame(index=[1,2,3,4,5],columns=["Observed","Expected"])
        both["Observed"] = real
        both["Expected"] = exp

        if MODEL == "linda":
            MODEL = "LINDA"
        elif MODEL == "linda-c":
            MODEL = "LINDA-C"
        elif MODEL == "deseq":
            MODEL = "DeSeq"
        elif MODEL == "deseq-c":
            MODEL = "DeSeq-C"
        elif MODEL == "ancom":
            MODEL = "ANCOM"
        elif MODEL == "ancom-c":
            MODEL = "ANCOM-C"
        elif MODEL == "aldex_we":
            MODEL = "ALDEx Welch"
        elif MODEL == "aldex_we-c":
            MODEL = "ALDEx Welch-C"
        elif MODEL == "aldex_wi":
            MODEL = "ALDEx Wilcoxon"
        elif MODEL == "aldex_wi-c":
            MODEL = "ALDEx Wilcoxon-C"
        elif MODEL == "ancom-bc2":
            MODEL = "ANCOM-BC2"
        elif MODEL == "ancom-bc2-c":
            MODEL = "ANCOM-BC2-C"
        elif MODEL == "lefse":
            MODEL = "LEfSe"
        elif MODEL =="mimic":
            MODEL = "miMic"
        elif MODEL == "mimic_relative":
            MODEL = "miMic relative"
        elif MODEL == "ancom-bc2_shuffled":
            MODEL = "ANCOM-BC2-(sh)"
    else:
        both = real
    if avg:
        both.plot(kind="bar", color=[color, "black"], rot=0, figsize=(5, 5), edgecolor="black")
    else:
        both.plot(kind="bar", color=[color, "black"], rot=0, figsize=(3.5, 3.5))
    plt.title(MODEL,fontsize=SIZE,fontweight="bold")
    plt.xlabel("No. of studies\nwhere significant", fontsize=SIZE)
    plt.ylabel("Percentage of\nsignificant specie", fontsize=SIZE)
    plt.xticks(fontsize=SIZE)
    plt.yticks(fontsize=SIZE)
    if avg:
        plt.legend(fontsize=17)
    else:
        plt.legend(fontsize=17,frameon=False)
    if avg:
        plt.ylim([0,80])
    plt.tight_layout()

    if avg and MODEL == "miMic":
        plt.savefig(f"rev_fig/consistency_of_model_between_datasets/IBD/real_exp_plots/avg_{MODEL}.png")
    else:
        plt.savefig(f"rev_fig/consistency_of_model_between_datasets/IBD/real_exp_plots/{MODEL}.png")
    plt.show()
    if not avg:
        return both
if __name__ == '__main__':
    mpl.rc('font', family='Times New Roman')
    LIST_MODELS = ["mimic_relative","mimic","linda", "linda-c", "deseq", "deseq-c", "ancom", "ancom-c", "aldex_we", "aldex_we-c", "aldex_wi",
                   "aldex_wi-c", "ancom-bc2", "ancom-bc2-c","ancom-bc2_shuffled", "lefse"]
    LIST_COLORS = ["mediumorchid","hotpink","lightskyblue","cornflowerblue","salmon","tomato","mediumaquamarine","seagreen","grey","dimgrey","lightslategray","slategray","green","darkgreen","green",'gold']

    list_both = list()
    for (MODEL,COLOR) in zip(LIST_MODELS,LIST_COLORS):
        real = pd.read_csv(f"rev_fig/consistency_of_model_between_datasets/IBD/real/{MODEL}.csv",index_col=0)
        exp = pd.read_csv(f"rev_fig/consistency_of_model_between_datasets/IBD/expected/{MODEL}.csv",index_col=0)
        both = plot_real_exp(real,exp,MODEL,COLOR)
        if MODEL == "mimic":
            plot_real_exp(both,both, "miMic", "hotpink", avg=True)

        if MODEL not in ["mimic_relative","mimic","ancom-bc2_shuffled"]:
            list_both.append(both)
    merged = pd.concat(list_both)
    avg_to_plot = merged.groupby(merged.index).mean()
    plot_real_exp(avg_to_plot, avg_to_plot, "Average over all models","white" ,avg=True)
    c=0
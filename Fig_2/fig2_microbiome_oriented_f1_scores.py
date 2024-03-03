import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
def calc_f1(df,data,model):
    TP = df.loc[data][f"TP-{model}"]
    FP = df.loc[data][f"FP-{model}"]
    FN = df.loc[data][f"FN-{model}"]
    F1 = 2*TP/(2*TP+FP+FN)
    return F1


def plot_f1(to_plot,data):
    SIZE = 20
    LIST_COLORS = ["hotpink", "lightskyblue", "cornflowerblue", "salmon", "tomato", "mediumaquamarine",
                   "seagreen", "grey", "dimgrey", "lightslategray", "slategray", "green", "darkgreen", 'gold']

    # Calculate the number of bars
    num_bars = len(to_plot)

    # Plot each bar with a unique color
    plt.figure(figsize=(4, 4))
    for i, (index, row) in enumerate(to_plot.iterrows()):

        color = LIST_COLORS[i % len(LIST_COLORS)]  # Ensure that colors repeat if there are more bars than colors
        plt.barh(index, row, color=color)

    plt.xlabel("F1 score",fontsize=SIZE)
    plt.xticks(fontsize=SIZE)
    plt.yticks(fontsize=SIZE)
    plt.title(data,fontsize=SIZE,fontweight="bold")
    plt.tight_layout()
    plt.savefig(f"ancom_bc2/{data}.png")
    plt.show()

def plot_f1_common(all,list_names):
    SIZE=25
    df1 = all[0]
    df2 = all[1]
    df3 = all[2]
    fig, axs = plt.subplots(1, 3, figsize=(15, 5), sharey=True)
    LIST_COLORS = ["hotpink", "lightskyblue", "cornflowerblue", "salmon", "tomato", "mediumaquamarine",
                   "seagreen", "grey", "dimgrey", "lightslategray", "slategray", "green", "darkgreen", 'gold']

    yticks_labels = LIST_MODELS
    axs[0].barh(df1.index, df1["F1"], color=LIST_COLORS)
    axs[1].barh(df2.index, df2["F1"], color=LIST_COLORS)
    axs[2].barh(df3.index, df2["F1"], color=LIST_COLORS)
    axs[0].set_xlabel("F1 score", fontsize=SIZE)
    axs[1].set_xlabel("F1 score", fontsize=SIZE)
    axs[2].set_xlabel("F1 score", fontsize=SIZE)

    # Add titles and labels
    axs[0].set_title(list_names[0], fontsize=30, fontweight="bold")
    axs[1].set_title(list_names[1], fontsize=30, fontweight="bold")
    axs[2].set_title(list_names[2], fontsize=30, fontweight="bold")
    axs[0].set_xticklabels([0.0,0.05,0.1,0.15,0.2,0.25], fontsize=SIZE)
    axs[1].set_xticklabels([0.0,0.1,0.2,0.3,0.4,0.5], fontsize=SIZE)
    axs[2].set_xticklabels([0.0,0.1,0.2,0.3,0.4,0.5], fontsize=SIZE)
    axs[0].set_yticklabels(yticks_labels, fontsize=20)

    plt.tight_layout()
    plt.savefig(f"ancom_bc2/simulations.png")
    plt.show()


if __name__ == '__main__':
    mpl.rc('font', family='Times New Roman')
    LIST_MODELS = ["miMic","LINDA", "LINDA-C", "DeSeq", "DeSeq-C", "ANCOM", "ANCOM-C", "ALDEx Welch",
                   "ALDEx Welch-C",
                   "ALDEx Wilcoxon", "ALDEx Wilcoxon-C", "ANCOM-BC2", "ANCOM-BC2-C","LEfSe"]
    df = pd.read_csv("ancom_bc2/all_sim_results.csv",index_col=0)
    #df = pd.read_csv("ancom_bc2/all_sim_results_common.csv",index_col=0)
    all_info = []
    for DATA in ["sim_ERP020401","sim_IBD","sim_PRJNA353587"]:
        to_plot = pd.DataFrame(index=LIST_MODELS,columns=["F1"])
        for MODEL in LIST_MODELS:
            to_plot["F1"][MODEL] = calc_f1(df,DATA,MODEL)
        all_info.append(to_plot)
        #plot_f1(to_plot, DATA)
    plot_f1_common(all_info,["Sim-ERP020401","Sim-IBD","Sim-PRJNA353587"])



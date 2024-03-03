import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.font_manager import FontProperties

if __name__ == '__main__':
    SIZE=25
    SIZE_T =30
    FS=(5,5)
    mpl.rc('font', family='Times New Roman')
    df1 = pd.read_csv("temp/outline_figures/Fig_2/tp_fp_00m_comp.csv",index_col=0)
    df2 = pd.read_csv("temp/outline_figures/Fig_2/fp_comp.csv",index_col=0)

##################################################################################################
    df1.plot(kind="barh", color=["hotpink", "lightpink", "deeppink"], rot=0, figsize=FS)
    plt.ylabel("Number of samples", fontsize=SIZE)
    plt.xlabel("TP", fontsize=SIZE)
    plt.title(r"Simulation - (0,0,$\mu$)", fontsize=SIZE_T, fontweight="bold")
    plt.xticks(fontsize=SIZE)
    plt.yticks(fontsize=SIZE)
    font = FontProperties()
    font.set_size(20)
    plt.legend(prop=font)
    plt.tight_layout()
    plt.savefig("temp/outline_figures/Fig_2/yoram/tp_0_0_m.png")
    plt.show()
    plt.clf()

#########################################################################################################
    df2.plot(kind="barh", color=["hotpink", "lightpink", "deeppink"], rot=0, figsize=FS,legend=False)
    plt.ylabel("Number of samples", fontsize=SIZE)
    plt.xlabel("FP", fontsize=SIZE)
    plt.title(r"Simulation - (0,0,$\mu$)", fontsize=SIZE_T, fontweight="bold")
    plt.xticks(fontsize=SIZE)
    plt.yticks(fontsize=SIZE)
    font = FontProperties()
    font.set_size(15)
    plt.tight_layout()
    plt.savefig("temp/outline_figures/Fig_2/yoram/fp_0_0_m.png")

    plt.show()

###########################################################################################################
    df = pd.read_csv("temp/outline_figures/Fig_2/sumulations_data/0_0_0_leaf_vs_mimic.csv",index_col=0)
    df.plot(kind="barh",color=["hotpink","lightpink","deeppink"],rot=0,figsize=FS,legend=False)
    plt.ylabel("Number of samples",fontsize=SIZE)
    plt.xlabel("FP",fontsize=SIZE)
    plt.title("Simulation - (0,0,0)",fontsize=SIZE_T,fontweight="bold")
    plt.xticks(fontsize=SIZE)
    plt.yticks(fontsize=SIZE)
    font = FontProperties()

    plt.tight_layout()
    plt.savefig("temp/outline_figures/Fig_2/yoram/fp_0_0_0.png")
    plt.show()

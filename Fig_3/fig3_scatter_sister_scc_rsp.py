import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
from scipy.stats import spearmanr

if __name__ == '__main__':
    WGS = True
    SIZE=15
    NAME="miMic-1"
    plt.figure(figsize=(4,4))
    mpl.rc('font', family='Times New Roman')
    df = pd.read_csv("temp/all_comparisons_16_7.csv",index_col=0)
    wgs = df.loc[["WGS-1","WGS-2","WGS-3","WGS-4","WGS-5","T2D","Obesity","Cirrhosis-POP"]]
    s16 = df.loc[["Cirrhosis","ERP020401","ERP021216","He","IBD","Jacob","MF","OK30","OK94","PRJNA353587","PRJNA419097","WB"]]
    scc3, p3 = spearmanr(df[NAME], df["SCC"])
    print(f"miMic-sis-first:\n SCC {scc3}\n p-value:{p3}")
    plt.scatter(wgs["SCC"], wgs[NAME],label=f"SCC:{round(scc3,3)} P-value:{round(p3,3)}",color="deeppink")
    plt.scatter(s16["SCC"], s16[NAME],color="hotpink")
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)
    plt.xlabel("Sister SCC",fontsize=15)
    plt.ylabel("RSP(1) score",fontsize=15)
    plt.legend()
    plt.tight_layout()
    plt.savefig("temp/outline_figures/Fig_3/sis_corr_scatter_1.png")
    plt.show()
    C=0
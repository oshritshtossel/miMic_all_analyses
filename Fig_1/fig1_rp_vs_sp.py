import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl

if __name__ == '__main__':
    SIZE = 15
    ALPHA =1
    MIMIC = True
    plt.figure(figsize=(4.5, 4.5))
    mpl.rc('font', family='Times New Roman')
    # Read results from folder
    df = pd.read_csv("Fig_1/Results/all_corrected_fp_tp.csv",index_col=0)

    # DeSeq
    plt.scatter(y=df["DeSeq-TP"]["Cirrhosis"],x=df["DeSeq-FP"]["Cirrhosis"],color="salmon",marker='X',alpha=ALPHA)
    plt.scatter(y=df["DeSeq-TP"]["Cirrhosis-POP"],x=df["DeSeq-FP"]["Cirrhosis-POP"],color="salmon",marker='v',alpha=ALPHA)
    plt.scatter(y=df["DeSeq-TP"]["ERP020401"],x=df["DeSeq-FP"]["ERP020401"],color="salmon",marker='+',alpha=ALPHA)
    plt.scatter(y=df["DeSeq-TP"]["ERP021216"],x=df["DeSeq-FP"]["ERP021216"],color="salmon",marker='D',alpha=ALPHA)
    plt.scatter(y=df["DeSeq-TP"]["He"],x=df["DeSeq-FP"]["He"],color="salmon",marker='d',alpha=ALPHA)
    plt.scatter(y=df["DeSeq-TP"]["IBD"],x=df["DeSeq-FP"]["IBD"],color="salmon",marker='_',alpha=ALPHA)
    plt.scatter(y=df["DeSeq-TP"]["Jacob"],x=df["DeSeq-FP"]["Jacob"],color="salmon",marker='|',alpha=ALPHA)
    plt.scatter(y=df["DeSeq-TP"]["MF"],x=df["DeSeq-FP"]["MF"],color="salmon",marker='x',alpha=ALPHA)
    plt.scatter(y=df["DeSeq-TP"]["Obesity"],x=df["DeSeq-FP"]["Obesity"],color="salmon",marker='^',alpha=ALPHA)
    plt.scatter(y=df["DeSeq-TP"]["OK30"],x=df["DeSeq-FP"]["OK30"],color="salmon",marker='*',alpha=ALPHA)
    plt.scatter(y=df["DeSeq-TP"]["OK94"],x=df["DeSeq-FP"]["OK94"],color="salmon",marker='.',alpha=ALPHA)
    plt.scatter(y=df["DeSeq-TP"]["PRJNA353587"],x=df["DeSeq-FP"]["PRJNA353587"],color="salmon",marker='s',alpha=ALPHA)
    plt.scatter(y=df["DeSeq-TP"]["PRJNA419097"],x=df["DeSeq-FP"]["PRJNA419097"],color="salmon",marker='p',alpha=ALPHA)
    plt.scatter(y=df["DeSeq-TP"]["T2D"],x=df["DeSeq-FP"]["T2D"],color="salmon",marker='8',alpha=ALPHA)
    plt.scatter(y=df["DeSeq-TP"]["WB"],x=df["DeSeq-FP"]["WB"],color="salmon",marker='h',alpha=ALPHA)
    plt.scatter(y=df["DeSeq-TP"]["WGS-1"],x=df["DeSeq-FP"]["WGS-1"],color="salmon",marker='<',alpha=ALPHA)
    plt.scatter(y=df["DeSeq-TP"]["WGS-2"],x=df["DeSeq-FP"]["WGS-2"],color="salmon",marker='>',alpha=ALPHA)
    plt.scatter(y=df["DeSeq-TP"]["WGS-3"],x=df["DeSeq-FP"]["WGS-3"],color="salmon",marker='4',alpha=ALPHA)
    plt.scatter(y=df["DeSeq-TP"]["WGS-4"],x=df["DeSeq-FP"]["WGS-4"],color="salmon",marker='3',alpha=ALPHA)
    plt.scatter(y=df["DeSeq-TP"]["WGS-5"],x=df["DeSeq-FP"]["WGS-5"],color="salmon",marker='2',alpha=ALPHA)

    # DeSeq-C
    plt.scatter(y=df["DeSeq-C-TP"]["Cirrhosis"], x=df["DeSeq-C-FP"]["Cirrhosis"], color="tomato", marker='X',alpha=ALPHA)
    plt.scatter(y=df["DeSeq-C-TP"]["Cirrhosis-POP"], x=df["DeSeq-C-FP"]["Cirrhosis-POP"], color="tomato", marker='v',alpha=ALPHA)
    plt.scatter(y=df["DeSeq-C-TP"]["ERP020401"], x=df["DeSeq-C-FP"]["ERP020401"], color="tomato", marker='+',alpha=ALPHA)
    plt.scatter(y=df["DeSeq-C-TP"]["ERP021216"], x=df["DeSeq-C-FP"]["ERP021216"], color="tomato", marker='D',alpha=ALPHA)
    plt.scatter(y=df["DeSeq-C-TP"]["He"], x=df["DeSeq-C-FP"]["He"], color="tomato", marker='d',alpha=ALPHA)
    plt.scatter(y=df["DeSeq-C-TP"]["IBD"], x=df["DeSeq-C-FP"]["IBD"], color="tomato", marker='_',alpha=ALPHA)
    plt.scatter(y=df["DeSeq-C-TP"]["Jacob"], x=df["DeSeq-C-FP"]["Jacob"], color="tomato", marker='|',alpha=ALPHA)
    plt.scatter(y=df["DeSeq-C-TP"]["MF"], x=df["DeSeq-C-FP"]["MF"], color="tomato", marker='x',alpha=ALPHA)
    plt.scatter(y=df["DeSeq-C-TP"]["Obesity"], x=df["DeSeq-C-FP"]["Obesity"], color="tomato", marker='^',alpha=ALPHA)
    plt.scatter(y=df["DeSeq-C-TP"]["OK30"], x=df["DeSeq-C-FP"]["OK30"], color="tomato", marker='*',alpha=ALPHA)
    plt.scatter(y=df["DeSeq-C-TP"]["OK94"], x=df["DeSeq-C-FP"]["OK94"], color="tomato", marker='.',alpha=ALPHA)
    plt.scatter(y=df["DeSeq-C-TP"]["PRJNA353587"], x=df["DeSeq-C-FP"]["PRJNA353587"], color="tomato", marker='s',alpha=ALPHA)
    plt.scatter(y=df["DeSeq-C-TP"]["PRJNA419097"], x=df["DeSeq-C-FP"]["PRJNA419097"], color="tomato", marker='p',alpha=ALPHA)
    plt.scatter(y=df["DeSeq-C-TP"]["T2D"], x=df["DeSeq-C-FP"]["T2D"], color="tomato", marker='8',alpha=ALPHA)
    plt.scatter(y=df["DeSeq-C-TP"]["WB"], x=df["DeSeq-C-FP"]["WB"], color="tomato", marker='h',alpha=ALPHA)
    plt.scatter(y=df["DeSeq-C-TP"]["WGS-1"], x=df["DeSeq-C-FP"]["WGS-1"], color="tomato", marker='<',alpha=ALPHA)
    plt.scatter(y=df["DeSeq-C-TP"]["WGS-2"], x=df["DeSeq-C-FP"]["WGS-2"], color="tomato", marker='>',alpha=ALPHA)
    plt.scatter(y=df["DeSeq-C-TP"]["WGS-3"], x=df["DeSeq-C-FP"]["WGS-3"], color="tomato", marker='4',alpha=ALPHA)
    plt.scatter(y=df["DeSeq-C-TP"]["WGS-4"], x=df["DeSeq-C-FP"]["WGS-4"], color="tomato", marker='3',alpha=ALPHA)
    plt.scatter(y=df["DeSeq-C-TP"]["WGS-5"], x=df["DeSeq-C-FP"]["WGS-5"], color="tomato", marker='2',alpha=ALPHA)

    # ANCOM
    plt.scatter(y=df["ANCOM-TP"]["Cirrhosis"], x=df["ANCOM-FP"]["Cirrhosis"], color="mediumaquamarine", marker='X',alpha=ALPHA)
    plt.scatter(y=df["ANCOM-TP"]["Cirrhosis-POP"], x=df["ANCOM-FP"]["Cirrhosis-POP"], color="mediumaquamarine", marker='v',alpha=ALPHA)
    plt.scatter(y=df["ANCOM-TP"]["ERP021216"], x=df["ANCOM-FP"]["ERP021216"], color="mediumaquamarine", marker='D',alpha=ALPHA)
    plt.scatter(y=df["ANCOM-TP"]["He"], x=df["ANCOM-FP"]["He"], color="mediumaquamarine", marker='d',alpha=ALPHA)
    plt.scatter(y=df["ANCOM-TP"]["IBD"], x=df["ANCOM-FP"]["IBD"], color="mediumaquamarine", marker='_',alpha=ALPHA)
    plt.scatter(y=df["ANCOM-TP"]["Jacob"], x=df["ANCOM-FP"]["Jacob"], color="mediumaquamarine", marker='|',alpha=ALPHA)
    plt.scatter(y=df["ANCOM-TP"]["MF"], x=df["ANCOM-FP"]["MF"], color="mediumaquamarine", marker='x',alpha=ALPHA)
    plt.scatter(y=df["ANCOM-TP"]["Obesity"], x=df["ANCOM-FP"]["Obesity"], color="mediumaquamarine", marker='^',alpha=ALPHA)
    plt.scatter(y=df["ANCOM-TP"]["OK30"], x=df["ANCOM-FP"]["OK30"], color="mediumaquamarine", marker='*',alpha=ALPHA)
    plt.scatter(y=df["ANCOM-TP"]["OK94"], x=df["ANCOM-FP"]["OK94"], color="mediumaquamarine", marker='.',alpha=ALPHA)
    plt.scatter(y=df["ANCOM-TP"]["PRJNA353587"], x=df["ANCOM-FP"]["PRJNA353587"], color="mediumaquamarine", marker='s',alpha=ALPHA)
    plt.scatter(y=df["ANCOM-TP"]["PRJNA419097"], x=df["ANCOM-FP"]["PRJNA419097"], color="mediumaquamarine", marker='p',alpha=ALPHA)
    plt.scatter(y=df["ANCOM-TP"]["T2D"], x=df["ANCOM-FP"]["T2D"], color="mediumaquamarine", marker='8',alpha=ALPHA)
    plt.scatter(y=df["ANCOM-TP"]["WB"], x=df["ANCOM-FP"]["WB"], color="mediumaquamarine", marker='h',alpha=ALPHA)
    plt.scatter(y=df["ANCOM-TP"]["WGS-1"], x=df["ANCOM-FP"]["WGS-1"], color="mediumaquamarine", marker='<',alpha=ALPHA)
    plt.scatter(y=df["ANCOM-TP"]["WGS-2"], x=df["ANCOM-FP"]["WGS-2"], color="mediumaquamarine", marker='>',alpha=ALPHA)
    plt.scatter(y=df["ANCOM-TP"]["WGS-3"], x=df["ANCOM-FP"]["WGS-3"], color="mediumaquamarine", marker='4',alpha=ALPHA)
    plt.scatter(y=df["ANCOM-TP"]["WGS-4"], x=df["ANCOM-FP"]["WGS-4"], color="mediumaquamarine", marker='3',alpha=ALPHA)
    plt.scatter(y=df["ANCOM-TP"]["WGS-5"], x=df["ANCOM-FP"]["WGS-5"], color="mediumaquamarine", marker='2',alpha=ALPHA)

    # ANCOM-C
    plt.scatter(y=df["ANCOM-C-TP"]["Cirrhosis"], x=df["ANCOM-C-FP"]["Cirrhosis"], color="seagreen", marker='X',alpha=ALPHA)
    plt.scatter(y=df["ANCOM-C-TP"]["Cirrhosis-POP"], x=df["ANCOM-C-FP"]["Cirrhosis-POP"], color="seagreen", marker='v',alpha=ALPHA)
    plt.scatter(y=df["ANCOM-C-TP"]["ERP020401"], x=df["ANCOM-C-FP"]["ERP020401"], color="seagreen", marker='+',alpha=ALPHA)
    plt.scatter(y=df["ANCOM-C-TP"]["ERP021216"], x=df["ANCOM-C-FP"]["ERP021216"], color="seagreen", marker='D',alpha=ALPHA)
    plt.scatter(y=df["ANCOM-C-TP"]["He"], x=df["ANCOM-C-FP"]["He"], color="seagreen", marker='d',alpha=ALPHA)
    plt.scatter(y=df["ANCOM-C-TP"]["IBD"], x=df["ANCOM-C-FP"]["IBD"], color="seagreen", marker='_',alpha=ALPHA)
    plt.scatter(y=df["ANCOM-C-TP"]["Jacob"], x=df["ANCOM-C-FP"]["Jacob"], color="seagreen", marker='|',alpha=ALPHA)
    plt.scatter(y=df["ANCOM-C-TP"]["MF"], x=df["ANCOM-C-FP"]["MF"], color="seagreen", marker='x',alpha=ALPHA)
    plt.scatter(y=df["ANCOM-C-TP"]["Obesity"], x=df["ANCOM-C-FP"]["Obesity"], color="seagreen", marker='^',alpha=ALPHA)
    plt.scatter(y=df["ANCOM-C-TP"]["OK30"], x=df["ANCOM-C-FP"]["OK30"], color="seagreen", marker='*',alpha=ALPHA)
    plt.scatter(y=df["ANCOM-C-TP"]["OK94"], x=df["ANCOM-C-FP"]["OK94"], color="seagreen", marker='.',alpha=ALPHA)
    plt.scatter(y=df["ANCOM-C-TP"]["PRJNA353587"], x=df["ANCOM-C-FP"]["PRJNA353587"], color="seagreen", marker='s',alpha=ALPHA)
    plt.scatter(y=df["ANCOM-C-TP"]["PRJNA419097"], x=df["ANCOM-C-FP"]["PRJNA419097"], color="seagreen", marker='p',alpha=ALPHA)
    plt.scatter(y=df["ANCOM-C-TP"]["T2D"], x=df["ANCOM-C-FP"]["T2D"], color="seagreen", marker='8',alpha=ALPHA)
    plt.scatter(y=df["ANCOM-C-TP"]["WB"], x=df["ANCOM-C-FP"]["WB"], color="seagreen", marker='h',alpha=ALPHA)
    plt.scatter(y=df["ANCOM-C-TP"]["WGS-1"], x=df["ANCOM-C-FP"]["WGS-1"], color="seagreen", marker='<',alpha=ALPHA)
    plt.scatter(y=df["ANCOM-C-TP"]["WGS-2"], x=df["ANCOM-C-FP"]["WGS-2"], color="seagreen", marker='>',alpha=ALPHA)
    plt.scatter(y=df["ANCOM-C-TP"]["WGS-3"], x=df["ANCOM-C-FP"]["WGS-3"], color="seagreen", marker='4',alpha=ALPHA)
    plt.scatter(y=df["ANCOM-C-TP"]["WGS-4"], x=df["ANCOM-C-FP"]["WGS-4"], color="seagreen", marker='3',alpha=ALPHA)
    plt.scatter(y=df["ANCOM-C-TP"]["WGS-5"], x=df["ANCOM-C-FP"]["WGS-5"], color="seagreen", marker='2',alpha=ALPHA)

    # ANCOM-BC2
    plt.scatter(y=df["TP-ANCOM-BC"]["Cirrhosis"], x=df["FP-ANCOM-BC"]["Cirrhosis"], color="green", marker='X',
                alpha=ALPHA)
    plt.scatter(y=df["TP-ANCOM-BC"]["Cirrhosis-POP"], x=df["FP-ANCOM-BC"]["Cirrhosis-POP"], color="green", marker='v',
                alpha=ALPHA)
    plt.scatter(y=df["TP-ANCOM-BC"]["ERP020401"], x=df["FP-ANCOM-BC"]["ERP020401"], color="green", marker='+',
                alpha=ALPHA)
    plt.scatter(y=df["TP-ANCOM-BC"]["ERP021216"], x=df["FP-ANCOM-BC"]["ERP021216"], color="green", marker='D',
                alpha=ALPHA)
    plt.scatter(y=df["TP-ANCOM-BC"]["He"], x=df["FP-ANCOM-BC"]["He"], color="green", marker='d', alpha=ALPHA)
    plt.scatter(y=df["TP-ANCOM-BC"]["IBD"], x=df["FP-ANCOM-BC"]["IBD"], color="green", marker='_', alpha=ALPHA)
    plt.scatter(y=df["TP-ANCOM-BC"]["Jacob"], x=df["FP-ANCOM-BC"]["Jacob"], color="green", marker='|', alpha=ALPHA)
    plt.scatter(y=df["TP-ANCOM-BC"]["MF"], x=df["FP-ANCOM-BC"]["MF"], color="green", marker='x', alpha=ALPHA)
    plt.scatter(y=df["TP-ANCOM-BC"]["Obesity"], x=df["FP-ANCOM-BC"]["Obesity"], color="green", marker='^', alpha=ALPHA)
    plt.scatter(y=df["TP-ANCOM-BC"]["OK30"], x=df["FP-ANCOM-BC"]["OK30"], color="green", marker='*', alpha=ALPHA)
    plt.scatter(y=df["TP-ANCOM-BC"]["OK94"], x=df["FP-ANCOM-BC"]["OK94"], color="green", marker='.', alpha=ALPHA)
    plt.scatter(y=df["TP-ANCOM-BC"]["PRJNA353587"], x=df["FP-ANCOM-BC"]["PRJNA353587"], color="green", marker='s',
                alpha=ALPHA)
    plt.scatter(y=df["TP-ANCOM-BC"]["PRJNA419097"], x=df["FP-ANCOM-BC"]["PRJNA419097"], color="green", marker='p',
                alpha=ALPHA)
    plt.scatter(y=df["TP-ANCOM-BC"]["T2D"], x=df["FP-ANCOM-BC"]["T2D"], color="green", marker='8', alpha=ALPHA)
    plt.scatter(y=df["TP-ANCOM-BC"]["WB"], x=df["FP-ANCOM-BC"]["WB"], color="green", marker='h', alpha=ALPHA)
    plt.scatter(y=df["TP-ANCOM-BC"]["WGS-1"], x=df["FP-ANCOM-BC"]["WGS-1"], color="green", marker='<', alpha=ALPHA)
    plt.scatter(y=df["TP-ANCOM-BC"]["WGS-2"], x=df["FP-ANCOM-BC"]["WGS-2"], color="green", marker='>', alpha=ALPHA)
    plt.scatter(y=df["TP-ANCOM-BC"]["WGS-3"], x=df["FP-ANCOM-BC"]["WGS-3"], color="green", marker='4', alpha=ALPHA)
    plt.scatter(y=df["TP-ANCOM-BC"]["WGS-4"], x=df["FP-ANCOM-BC"]["WGS-4"], color="green", marker='3', alpha=ALPHA)
    plt.scatter(y=df["TP-ANCOM-BC"]["WGS-5"], x=df["FP-ANCOM-BC"]["WGS-5"], color="green", marker='2', alpha=ALPHA)


    # ANCOM-BC2-C
    plt.scatter(y=df["TP-ANCOM-BC-C"]["Cirrhosis"], x=df["FP-ANCOM-BC-C"]["Cirrhosis"], color="darkgreen", marker='X',
                alpha=ALPHA)
    plt.scatter(y=df["TP-ANCOM-BC-C"]["Cirrhosis-POP"], x=df["FP-ANCOM-BC-C"]["Cirrhosis-POP"], color="darkgreen",
                marker='v',
                alpha=ALPHA)
    plt.scatter(y=df["TP-ANCOM-BC-C"]["ERP020401"], x=df["FP-ANCOM-BC-C"]["ERP020401"], color="darkgreen", marker='+',
                alpha=ALPHA)
    plt.scatter(y=df["TP-ANCOM-BC-C"]["ERP021216"], x=df["FP-ANCOM-BC-C"]["ERP021216"], color="darkgreen", marker='D',
                alpha=ALPHA)
    plt.scatter(y=df["TP-ANCOM-BC-C"]["He"], x=df["FP-ANCOM-BC-C"]["He"], color="darkgreen", marker='d', alpha=ALPHA)
    plt.scatter(y=df["TP-ANCOM-BC-C"]["IBD"], x=df["FP-ANCOM-BC-C"]["IBD"], color="darkgreen", marker='_', alpha=ALPHA)
    plt.scatter(y=df["TP-ANCOM-BC-C"]["Jacob"], x=df["FP-ANCOM-BC-C"]["Jacob"], color="darkgreen", marker='|', alpha=ALPHA)
    plt.scatter(y=df["TP-ANCOM-BC-C"]["MF"], x=df["FP-ANCOM-BC-C"]["MF"], color="darkgreen", marker='x', alpha=ALPHA)
    plt.scatter(y=df["TP-ANCOM-BC-C"]["Obesity"], x=df["FP-ANCOM-BC-C"]["Obesity"], color="darkgreen", marker='^', alpha=ALPHA)
    plt.scatter(y=df["TP-ANCOM-BC-C"]["OK30"], x=df["FP-ANCOM-BC-C"]["OK30"], color="darkgreen", marker='*', alpha=ALPHA)
    plt.scatter(y=df["TP-ANCOM-BC-C"]["OK94"], x=df["FP-ANCOM-BC-C"]["OK94"], color="darkgreen", marker='.', alpha=ALPHA)
    plt.scatter(y=df["TP-ANCOM-BC-C"]["PRJNA353587"], x=df["FP-ANCOM-BC-C"]["PRJNA353587"], color="darkgreen", marker='s',
                alpha=ALPHA)
    plt.scatter(y=df["TP-ANCOM-BC-C"]["PRJNA419097"], x=df["FP-ANCOM-BC-C"]["PRJNA419097"], color="darkgreen", marker='p',
                alpha=ALPHA)
    plt.scatter(y=df["TP-ANCOM-BC-C"]["T2D"], x=df["FP-ANCOM-BC-C"]["T2D"], color="darkgreen", marker='8', alpha=ALPHA)
    plt.scatter(y=df["TP-ANCOM-BC-C"]["WB"], x=df["FP-ANCOM-BC-C"]["WB"], color="darkgreen", marker='h', alpha=ALPHA)
    plt.scatter(y=df["TP-ANCOM-BC-C"]["WGS-1"], x=df["FP-ANCOM-BC-C"]["WGS-1"], color="darkgreen", marker='<', alpha=ALPHA)
    plt.scatter(y=df["TP-ANCOM-BC-C"]["WGS-2"], x=df["FP-ANCOM-BC-C"]["WGS-2"], color="darkgreen", marker='>', alpha=ALPHA)
    plt.scatter(y=df["TP-ANCOM-BC-C"]["WGS-3"], x=df["FP-ANCOM-BC-C"]["WGS-3"], color="darkgreen", marker='4', alpha=ALPHA)
    plt.scatter(y=df["TP-ANCOM-BC-C"]["WGS-4"], x=df["FP-ANCOM-BC-C"]["WGS-4"], color="darkgreen", marker='3', alpha=ALPHA)
    plt.scatter(y=df["TP-ANCOM-BC-C"]["WGS-5"], x=df["FP-ANCOM-BC-C"]["WGS-5"], color="darkgreen", marker='2', alpha=ALPHA)

    # LINDA
    plt.scatter(y=df["LINDA-TP"]["Cirrhosis"], x=df["LINDA-FP"]["Cirrhosis"], color="lightskyblue", marker='X',alpha=ALPHA)
    plt.scatter(y=df["LINDA-TP"]["Cirrhosis-POP"], x=df["LINDA-FP"]["Cirrhosis-POP"], color="lightskyblue", marker='v',alpha=ALPHA)
    plt.scatter(y=df["LINDA-TP"]["ERP020401"], x=df["LINDA-FP"]["ERP020401"], color="lightskyblue", marker='+',alpha=ALPHA)
    plt.scatter(y=df["LINDA-TP"]["ERP021216"], x=df["LINDA-FP"]["ERP021216"], color="lightskyblue", marker='D',alpha=ALPHA)
    plt.scatter(y=df["LINDA-TP"]["He"], x=df["LINDA-FP"]["He"], color="lightskyblue", marker='d',alpha=ALPHA)
    plt.scatter(y=df["LINDA-TP"]["IBD"], x=df["LINDA-FP"]["IBD"], color="lightskyblue", marker='_',alpha=ALPHA)
    plt.scatter(y=df["LINDA-TP"]["Jacob"], x=df["LINDA-FP"]["Jacob"], color="lightskyblue", marker='|',alpha=ALPHA)
    plt.scatter(y=df["LINDA-TP"]["MF"], x=df["LINDA-FP"]["MF"], color="lightskyblue", marker='x',alpha=ALPHA)
    plt.scatter(y=df["LINDA-TP"]["Obesity"], x=df["LINDA-FP"]["Obesity"], color="lightskyblue", marker='^',alpha=ALPHA)
    plt.scatter(y=df["LINDA-TP"]["OK30"], x=df["LINDA-FP"]["OK30"], color="lightskyblue", marker='*',alpha=ALPHA)
    plt.scatter(y=df["LINDA-TP"]["OK94"], x=df["LINDA-FP"]["OK94"], color="lightskyblue", marker='.',alpha=ALPHA)
    plt.scatter(y=df["LINDA-TP"]["PRJNA353587"], x=df["LINDA-FP"]["PRJNA353587"], color="lightskyblue", marker='s',alpha=ALPHA)
    plt.scatter(y=df["LINDA-TP"]["PRJNA419097"], x=df["LINDA-FP"]["PRJNA419097"], color="lightskyblue", marker='p',alpha=ALPHA)
    plt.scatter(y=df["LINDA-TP"]["T2D"], x=df["LINDA-FP"]["T2D"], color="lightskyblue", marker='8',alpha=ALPHA)
    plt.scatter(y=df["LINDA-TP"]["WB"], x=df["LINDA-FP"]["WB"], color="lightskyblue", marker='h',alpha=ALPHA)
    plt.scatter(y=df["LINDA-TP"]["WGS-1"], x=df["LINDA-FP"]["WGS-1"], color="lightskyblue", marker='<',alpha=ALPHA)
    plt.scatter(y=df["LINDA-TP"]["WGS-2"], x=df["LINDA-FP"]["WGS-2"], color="lightskyblue", marker='>',alpha=ALPHA)
    plt.scatter(y=df["LINDA-TP"]["WGS-3"], x=df["LINDA-FP"]["WGS-3"], color="lightskyblue", marker='4',alpha=ALPHA)
    plt.scatter(y=df["LINDA-TP"]["WGS-4"], x=df["LINDA-FP"]["WGS-4"], color="lightskyblue", marker='3',alpha=ALPHA)
    plt.scatter(y=df["LINDA-TP"]["WGS-5"], x=df["LINDA-FP"]["WGS-5"], color="lightskyblue", marker='2',alpha=ALPHA)

    # LINDA-C
    plt.scatter(y=df["LINDA-C-TP"]["Cirrhosis"], x=df["LINDA-C-FP"]["Cirrhosis"], color="cornflowerblue", marker='X',alpha=ALPHA)
    plt.scatter(y=df["LINDA-C-TP"]["Cirrhosis-POP"], x=df["LINDA-C-FP"]["Cirrhosis-POP"], color="cornflowerblue", marker='v',alpha=ALPHA)
    plt.scatter(y=df["LINDA-C-TP"]["ERP020401"], x=df["LINDA-C-FP"]["ERP020401"], color="cornflowerblue", marker='+',alpha=ALPHA)
    plt.scatter(y=df["LINDA-C-TP"]["ERP021216"], x=df["LINDA-C-FP"]["ERP021216"], color="cornflowerblue", marker='D',alpha=ALPHA)
    plt.scatter(y=df["LINDA-C-TP"]["He"], x=df["LINDA-C-FP"]["He"], color="cornflowerblue", marker='d',alpha=ALPHA)
    plt.scatter(y=df["LINDA-C-TP"]["IBD"], x=df["LINDA-C-FP"]["IBD"], color="cornflowerblue", marker='_',alpha=ALPHA)
    plt.scatter(y=df["LINDA-C-TP"]["Jacob"], x=df["LINDA-C-FP"]["Jacob"], color="cornflowerblue", marker='|',alpha=ALPHA)
    plt.scatter(y=df["LINDA-C-TP"]["MF"], x=df["LINDA-C-FP"]["MF"], color="cornflowerblue", marker='x',alpha=ALPHA)
    plt.scatter(y=df["LINDA-C-TP"]["Obesity"], x=df["LINDA-C-FP"]["Obesity"], color="cornflowerblue", marker='^',alpha=ALPHA)
    plt.scatter(y=df["LINDA-C-TP"]["OK30"], x=df["LINDA-C-FP"]["OK30"], color="cornflowerblue", marker='*',alpha=ALPHA)
    plt.scatter(y=df["LINDA-C-TP"]["OK94"], x=df["LINDA-C-FP"]["OK94"], color="cornflowerblue", marker='.',alpha=ALPHA)
    plt.scatter(y=df["LINDA-C-TP"]["PRJNA353587"], x=df["LINDA-C-FP"]["PRJNA353587"], color="cornflowerblue", marker='s',alpha=ALPHA)
    plt.scatter(y=df["LINDA-C-TP"]["PRJNA419097"], x=df["LINDA-C-FP"]["PRJNA419097"], color="cornflowerblue", marker='p',alpha=ALPHA)
    plt.scatter(y=df["LINDA-C-TP"]["T2D"], x=df["LINDA-C-FP"]["T2D"], color="cornflowerblue", marker='8',alpha=ALPHA)
    plt.scatter(y=df["LINDA-C-TP"]["WB"], x=df["LINDA-C-FP"]["WB"], color="cornflowerblue", marker='h',alpha=ALPHA)
    plt.scatter(y=df["LINDA-C-TP"]["WGS-1"], x=df["LINDA-C-FP"]["WGS-1"], color="cornflowerblue", marker='<',alpha=ALPHA)
    plt.scatter(y=df["LINDA-C-TP"]["WGS-2"], x=df["LINDA-C-FP"]["WGS-2"], color="cornflowerblue", marker='>',alpha=ALPHA)
    plt.scatter(y=df["LINDA-C-TP"]["WGS-3"], x=df["LINDA-C-FP"]["WGS-3"], color="cornflowerblue", marker='4',alpha=ALPHA)
    plt.scatter(y=df["LINDA-C-TP"]["WGS-4"], x=df["LINDA-C-FP"]["WGS-4"], color="cornflowerblue", marker='3',alpha=ALPHA)
    plt.scatter(y=df["LINDA-C-TP"]["WGS-5"], x=df["LINDA-C-FP"]["WGS-5"], color="cornflowerblue", marker='2',alpha=ALPHA)

    # ALDEx Welch
    plt.scatter(y=df["TP-ALDEx Welch"]["Cirrhosis"], x=df["FP-ALDEx Welch"]["Cirrhosis"], color="grey", marker='X',
                alpha=ALPHA)
    plt.scatter(y=df["TP-ALDEx Welch"]["Cirrhosis-POP"], x=df["FP-ALDEx Welch"]["Cirrhosis-POP"], color="grey",
                marker='v', alpha=ALPHA)
    plt.scatter(y=df["TP-ALDEx Welch"]["ERP020401"], x=df["FP-ALDEx Welch"]["ERP020401"], color="grey", marker='+',
                alpha=ALPHA)
    plt.scatter(y=df["TP-ALDEx Welch"]["ERP021216"], x=df["FP-ALDEx Welch"]["ERP021216"], color="grey", marker='D',
                alpha=ALPHA)
    plt.scatter(y=df["TP-ALDEx Welch"]["He"], x=df["FP-ALDEx Welch"]["He"], color="grey", marker='d', alpha=ALPHA)
    plt.scatter(y=df["TP-ALDEx Welch"]["IBD"], x=df["FP-ALDEx Welch"]["IBD"], color="grey", marker='_', alpha=ALPHA)
    plt.scatter(y=df["TP-ALDEx Welch"]["Jacob"], x=df["FP-ALDEx Welch"]["Jacob"], color="grey", marker='|',
                alpha=ALPHA)
    plt.scatter(y=df["TP-ALDEx Welch"]["MF"], x=df["FP-ALDEx Welch"]["MF"], color="grey", marker='x', alpha=ALPHA)
    plt.scatter(y=df["TP-ALDEx Welch"]["Obesity"], x=df["FP-ALDEx Welch"]["Obesity"], color="grey", marker='^',
                alpha=ALPHA)
    plt.scatter(y=df["TP-ALDEx Welch"]["OK30"], x=df["FP-ALDEx Welch"]["OK30"], color="grey", marker='*', alpha=ALPHA)
    plt.scatter(y=df["TP-ALDEx Welch"]["OK94"], x=df["FP-ALDEx Welch"]["OK94"], color="grey", marker='.', alpha=ALPHA)
    plt.scatter(y=df["TP-ALDEx Welch"]["PRJNA353587"], x=df["FP-ALDEx Welch"]["PRJNA353587"], color="grey",
                marker='s', alpha=ALPHA)
    plt.scatter(y=df["TP-ALDEx Welch"]["PRJNA419097"], x=df["FP-ALDEx Welch"]["PRJNA419097"], color="grey",
                marker='p', alpha=ALPHA)
    plt.scatter(y=df["TP-ALDEx Welch"]["T2D"], x=df["FP-ALDEx Welch"]["T2D"], color="grey", marker='8', alpha=ALPHA)
    plt.scatter(y=df["TP-ALDEx Welch"]["WB"], x=df["FP-ALDEx Welch"]["WB"], color="grey", marker='h', alpha=ALPHA)
    plt.scatter(y=df["TP-ALDEx Welch"]["WGS-1"], x=df["FP-ALDEx Welch"]["WGS-1"], color="grey", marker='<',
                alpha=ALPHA)
    plt.scatter(y=df["TP-ALDEx Welch"]["WGS-2"], x=df["FP-ALDEx Welch"]["WGS-2"], color="grey", marker='>',
                alpha=ALPHA)
    plt.scatter(y=df["TP-ALDEx Welch"]["WGS-3"], x=df["FP-ALDEx Welch"]["WGS-3"], color="grey", marker='4',
                alpha=ALPHA)
    plt.scatter(y=df["TP-ALDEx Welch"]["WGS-4"], x=df["FP-ALDEx Welch"]["WGS-4"], color="grey", marker='3',
                alpha=ALPHA)
    plt.scatter(y=df["TP-ALDEx Welch"]["WGS-5"], x=df["FP-ALDEx Welch"]["WGS-5"], color="grey", marker='2',
                alpha=ALPHA)

    # ALDEx Welch-C
    plt.scatter(y=df["TP-ALDEx-c Welch"]["Cirrhosis"], x=df["FP-ALDEx-c Welch"]["Cirrhosis"], color="dimgrey", marker='X',
                alpha=ALPHA)
    plt.scatter(y=df["TP-ALDEx-c Welch"]["Cirrhosis-POP"], x=df["FP-ALDEx-c Welch"]["Cirrhosis-POP"], color="dimgrey",
                marker='v', alpha=ALPHA)
    plt.scatter(y=df["TP-ALDEx-c Welch"]["ERP020401"], x=df["FP-ALDEx-c Welch"]["ERP020401"], color="dimgrey", marker='+',
                alpha=ALPHA)
    plt.scatter(y=df["TP-ALDEx-c Welch"]["ERP021216"], x=df["FP-ALDEx-c Welch"]["ERP021216"], color="dimgrey", marker='D',
                alpha=ALPHA)
    plt.scatter(y=df["TP-ALDEx-c Welch"]["He"], x=df["FP-ALDEx-c Welch"]["He"], color="dimgrey", marker='d', alpha=ALPHA)
    plt.scatter(y=df["TP-ALDEx-c Welch"]["IBD"], x=df["FP-ALDEx-c Welch"]["IBD"], color="dimgrey", marker='_', alpha=ALPHA)
    plt.scatter(y=df["TP-ALDEx-c Welch"]["Jacob"], x=df["FP-ALDEx-c Welch"]["Jacob"], color="dimgrey", marker='|',
                alpha=ALPHA)
    plt.scatter(y=df["TP-ALDEx-c Welch"]["MF"], x=df["FP-ALDEx-c Welch"]["MF"], color="dimgrey", marker='x', alpha=ALPHA)
    plt.scatter(y=df["TP-ALDEx-c Welch"]["Obesity"], x=df["FP-ALDEx-c Welch"]["Obesity"], color="dimgrey", marker='^',
                alpha=ALPHA)
    plt.scatter(y=df["TP-ALDEx-c Welch"]["OK30"], x=df["FP-ALDEx-c Welch"]["OK30"], color="dimgrey", marker='*', alpha=ALPHA)
    plt.scatter(y=df["TP-ALDEx-c Welch"]["OK94"], x=df["FP-ALDEx-c Welch"]["OK94"], color="dimgrey", marker='.', alpha=ALPHA)
    plt.scatter(y=df["TP-ALDEx-c Welch"]["PRJNA353587"], x=df["FP-ALDEx-c Welch"]["PRJNA353587"], color="dimgrey",
                marker='s', alpha=ALPHA)
    plt.scatter(y=df["TP-ALDEx-c Welch"]["PRJNA419097"], x=df["FP-ALDEx-c Welch"]["PRJNA419097"], color="dimgrey",
                marker='p', alpha=ALPHA)
    plt.scatter(y=df["TP-ALDEx-c Welch"]["T2D"], x=df["FP-ALDEx-c Welch"]["T2D"], color="dimgrey", marker='8', alpha=ALPHA)
    plt.scatter(y=df["TP-ALDEx-c Welch"]["WB"], x=df["FP-ALDEx-c Welch"]["WB"], color="dimgrey", marker='h', alpha=ALPHA)
    plt.scatter(y=df["TP-ALDEx-c Welch"]["WGS-1"], x=df["FP-ALDEx-c Welch"]["WGS-1"], color="dimgrey", marker='<',
                alpha=ALPHA)
    plt.scatter(y=df["TP-ALDEx-c Welch"]["WGS-2"], x=df["FP-ALDEx-c Welch"]["WGS-2"], color="dimgrey", marker='>',
                alpha=ALPHA)
    plt.scatter(y=df["TP-ALDEx-c Welch"]["WGS-3"], x=df["FP-ALDEx-c Welch"]["WGS-3"], color="dimgrey", marker='4',
                alpha=ALPHA)
    plt.scatter(y=df["TP-ALDEx-c Welch"]["WGS-4"], x=df["FP-ALDEx-c Welch"]["WGS-4"], color="dimgrey", marker='3',
                alpha=ALPHA)
    plt.scatter(y=df["TP-ALDEx-c Welch"]["WGS-5"], x=df["FP-ALDEx-c Welch"]["WGS-5"], color="dimgrey", marker='2',
                alpha=ALPHA)

    # ALDEx Wilcoxon
    plt.scatter(y=df["TP-ALDEx Wilcoxon"]["Cirrhosis"], x=df["FP-ALDEx Wilcoxon"]["Cirrhosis"], color="lightslategray", marker='X',
                alpha=ALPHA)
    plt.scatter(y=df["TP-ALDEx Wilcoxon"]["Cirrhosis-POP"], x=df["FP-ALDEx Wilcoxon"]["Cirrhosis-POP"], color="lightslategray", marker='v',
                alpha=ALPHA)
    plt.scatter(y=df["TP-ALDEx Wilcoxon"]["ERP020401"], x=df["FP-ALDEx Wilcoxon"]["ERP020401"], color="lightslategray", marker='+',
                alpha=ALPHA)
    plt.scatter(y=df["TP-ALDEx Wilcoxon"]["ERP021216"], x=df["FP-ALDEx Wilcoxon"]["ERP021216"], color="lightslategray", marker='D',
                alpha=ALPHA)
    plt.scatter(y=df["TP-ALDEx Wilcoxon"]["He"], x=df["FP-ALDEx Wilcoxon"]["He"], color="lightslategray", marker='d', alpha=ALPHA)
    plt.scatter(y=df["TP-ALDEx Wilcoxon"]["IBD"], x=df["FP-ALDEx Wilcoxon"]["IBD"], color="lightslategray", marker='_', alpha=ALPHA)
    plt.scatter(y=df["TP-ALDEx Wilcoxon"]["Jacob"], x=df["FP-ALDEx Wilcoxon"]["Jacob"], color="lightslategray", marker='|', alpha=ALPHA)
    plt.scatter(y=df["TP-ALDEx Wilcoxon"]["MF"], x=df["FP-ALDEx Wilcoxon"]["MF"], color="lightslategray", marker='x', alpha=ALPHA)
    plt.scatter(y=df["TP-ALDEx Wilcoxon"]["Obesity"], x=df["FP-ALDEx Wilcoxon"]["Obesity"], color="lightslategray", marker='^', alpha=ALPHA)
    plt.scatter(y=df["TP-ALDEx Wilcoxon"]["OK30"], x=df["FP-ALDEx Wilcoxon"]["OK30"], color="lightslategray", marker='*', alpha=ALPHA)
    plt.scatter(y=df["TP-ALDEx Wilcoxon"]["OK94"], x=df["FP-ALDEx Wilcoxon"]["OK94"], color="lightslategray", marker='.', alpha=ALPHA)
    plt.scatter(y=df["TP-ALDEx Wilcoxon"]["PRJNA353587"], x=df["FP-ALDEx Wilcoxon"]["PRJNA353587"], color="lightslategray", marker='s',
                alpha=ALPHA)
    plt.scatter(y=df["TP-ALDEx Wilcoxon"]["PRJNA419097"], x=df["FP-ALDEx Wilcoxon"]["PRJNA419097"], color="lightslategray", marker='p',
                alpha=ALPHA)
    plt.scatter(y=df["TP-ALDEx Wilcoxon"]["T2D"], x=df["FP-ALDEx Wilcoxon"]["T2D"], color="lightslategray", marker='8', alpha=ALPHA)
    plt.scatter(y=df["TP-ALDEx Wilcoxon"]["WB"], x=df["FP-ALDEx Wilcoxon"]["WB"], color="lightslategray", marker='h', alpha=ALPHA)
    plt.scatter(y=df["TP-ALDEx Wilcoxon"]["WGS-1"], x=df["FP-ALDEx Wilcoxon"]["WGS-1"], color="lightslategray", marker='<', alpha=ALPHA)
    plt.scatter(y=df["TP-ALDEx Wilcoxon"]["WGS-2"], x=df["FP-ALDEx Wilcoxon"]["WGS-2"], color="lightslategray", marker='>', alpha=ALPHA)
    plt.scatter(y=df["TP-ALDEx Wilcoxon"]["WGS-3"], x=df["FP-ALDEx Wilcoxon"]["WGS-3"], color="lightslategray", marker='4', alpha=ALPHA)
    plt.scatter(y=df["TP-ALDEx Wilcoxon"]["WGS-4"], x=df["FP-ALDEx Wilcoxon"]["WGS-4"], color="lightslategray", marker='3', alpha=ALPHA)
    plt.scatter(y=df["TP-ALDEx Wilcoxon"]["WGS-5"], x=df["FP-ALDEx Wilcoxon"]["WGS-5"], color="lightslategray", marker='2', alpha=ALPHA)

    # ALDEx Wilcoxon-C
    plt.scatter(y=df["TP-ALDEx-c Wilcoxon"]["Cirrhosis"], x=df["FP-ALDEx-c Wilcoxon"]["Cirrhosis"], color="slategray",
                marker='X',
                alpha=ALPHA)
    plt.scatter(y=df["TP-ALDEx-c Wilcoxon"]["Cirrhosis-POP"], x=df["FP-ALDEx-c Wilcoxon"]["Cirrhosis-POP"], color="slategray",
                marker='v',
                alpha=ALPHA)
    plt.scatter(y=df["TP-ALDEx-c Wilcoxon"]["ERP020401"], x=df["FP-ALDEx-c Wilcoxon"]["ERP020401"], color="slategray",
                marker='+',
                alpha=ALPHA)
    plt.scatter(y=df["TP-ALDEx-c Wilcoxon"]["ERP021216"], x=df["FP-ALDEx-c Wilcoxon"]["ERP021216"], color="slategray",
                marker='D',
                alpha=ALPHA)
    plt.scatter(y=df["TP-ALDEx-c Wilcoxon"]["He"], x=df["FP-ALDEx-c Wilcoxon"]["He"], color="slategray", marker='d',
                alpha=ALPHA)
    plt.scatter(y=df["TP-ALDEx-c Wilcoxon"]["IBD"], x=df["FP-ALDEx-c Wilcoxon"]["IBD"], color="slategray", marker='_',
                alpha=ALPHA)
    plt.scatter(y=df["TP-ALDEx-c Wilcoxon"]["Jacob"], x=df["FP-ALDEx-c Wilcoxon"]["Jacob"], color="slategray",
                marker='|', alpha=ALPHA)
    plt.scatter(y=df["TP-ALDEx-c Wilcoxon"]["MF"], x=df["FP-ALDEx-c Wilcoxon"]["MF"], color="slategray", marker='x',
                alpha=ALPHA)
    plt.scatter(y=df["TP-ALDEx-c Wilcoxon"]["Obesity"], x=df["FP-ALDEx-c Wilcoxon"]["Obesity"], color="slategray",
                marker='^', alpha=ALPHA)
    plt.scatter(y=df["TP-ALDEx-c Wilcoxon"]["OK30"], x=df["FP-ALDEx-c Wilcoxon"]["OK30"], color="slategray",
                marker='*', alpha=ALPHA)
    plt.scatter(y=df["TP-ALDEx-c Wilcoxon"]["OK94"], x=df["FP-ALDEx-c Wilcoxon"]["OK94"], color="slategray",
                marker='.', alpha=ALPHA)
    plt.scatter(y=df["TP-ALDEx-c Wilcoxon"]["PRJNA353587"], x=df["FP-ALDEx-c Wilcoxon"]["PRJNA353587"],
                color="slategray", marker='s',
                alpha=ALPHA)
    plt.scatter(y=df["TP-ALDEx-c Wilcoxon"]["PRJNA419097"], x=df["FP-ALDEx-c Wilcoxon"]["PRJNA419097"],
                color="slategray", marker='p',
                alpha=ALPHA)
    plt.scatter(y=df["TP-ALDEx-c Wilcoxon"]["T2D"], x=df["FP-ALDEx-c Wilcoxon"]["T2D"], color="slategray", marker='8',
                alpha=ALPHA)
    plt.scatter(y=df["TP-ALDEx-c Wilcoxon"]["WB"], x=df["FP-ALDEx-c Wilcoxon"]["WB"], color="slategray", marker='h',
                alpha=ALPHA)
    plt.scatter(y=df["TP-ALDEx-c Wilcoxon"]["WGS-1"], x=df["FP-ALDEx-c Wilcoxon"]["WGS-1"], color="slategray",
                marker='<', alpha=ALPHA)
    plt.scatter(y=df["TP-ALDEx-c Wilcoxon"]["WGS-2"], x=df["FP-ALDEx-c Wilcoxon"]["WGS-2"], color="slategray",
                marker='>', alpha=ALPHA)
    plt.scatter(y=df["TP-ALDEx-c Wilcoxon"]["WGS-3"], x=df["FP-ALDEx-c Wilcoxon"]["WGS-3"], color="slategray",
                marker='4', alpha=ALPHA)
    plt.scatter(y=df["TP-ALDEx-c Wilcoxon"]["WGS-4"], x=df["FP-ALDEx-c Wilcoxon"]["WGS-4"], color="slategray",
                marker='3', alpha=ALPHA)
    plt.scatter(y=df["TP-ALDEx-c Wilcoxon"]["WGS-5"], x=df["FP-ALDEx-c Wilcoxon"]["WGS-5"], color="slategray",
                marker='2', alpha=ALPHA)

    # LEfSe
    plt.scatter(y=df["LEfSe-TP"]["Cirrhosis"], x=df["LEfSe-FP"]["Cirrhosis"], color="gold", marker='X',
                alpha=ALPHA)
    plt.scatter(y=df["LEfSe-TP"]["Cirrhosis-POP"], x=df["LEfSe-FP"]["Cirrhosis-POP"], color="gold", marker='v',
                alpha=ALPHA)
    plt.scatter(y=df["LEfSe-TP"]["ERP020401"], x=df["LEfSe-FP"]["ERP020401"], color="gold", marker='+',
                alpha=ALPHA)
    plt.scatter(y=df["LEfSe-TP"]["ERP021216"], x=df["LEfSe-FP"]["ERP021216"], color="gold", marker='D',
                alpha=ALPHA)
    plt.scatter(y=df["LEfSe-TP"]["He"], x=df["LEfSe-FP"]["He"], color="gold", marker='d', alpha=ALPHA)
    plt.scatter(y=df["LEfSe-TP"]["IBD"], x=df["LEfSe-FP"]["IBD"], color="gold", marker='_', alpha=ALPHA)
    plt.scatter(y=df["LEfSe-TP"]["Jacob"], x=df["LEfSe-FP"]["Jacob"], color="gold", marker='|', alpha=ALPHA)
    plt.scatter(y=df["LEfSe-TP"]["MF"], x=df["LEfSe-FP"]["MF"], color="gold", marker='x', alpha=ALPHA)
    plt.scatter(y=df["LEfSe-TP"]["Obesity"], x=df["LEfSe-FP"]["Obesity"], color="gold", marker='^', alpha=ALPHA)
    plt.scatter(y=df["LEfSe-TP"]["OK30"], x=df["LEfSe-FP"]["OK30"], color="gold", marker='*', alpha=ALPHA)
    plt.scatter(y=df["LEfSe-TP"]["OK94"], x=df["LEfSe-FP"]["OK94"], color="gold", marker='.', alpha=ALPHA)
    plt.scatter(y=df["LEfSe-TP"]["PRJNA353587"], x=df["LEfSe-FP"]["PRJNA353587"], color="gold", marker='s',
                alpha=ALPHA)
    plt.scatter(y=df["LEfSe-TP"]["PRJNA419097"], x=df["LEfSe-FP"]["PRJNA419097"], color="gold", marker='p',
                alpha=ALPHA)
    plt.scatter(y=df["LEfSe-TP"]["T2D"], x=df["LEfSe-FP"]["T2D"], color="gold", marker='8', alpha=ALPHA)
    plt.scatter(y=df["LEfSe-TP"]["WB"], x=df["LEfSe-FP"]["WB"], color="gold", marker='h', alpha=ALPHA)
    plt.scatter(y=df["LEfSe-TP"]["WGS-1"], x=df["LEfSe-FP"]["WGS-1"], color="gold", marker='<', alpha=ALPHA)
    plt.scatter(y=df["LEfSe-TP"]["WGS-2"], x=df["LEfSe-FP"]["WGS-2"], color="gold", marker='>', alpha=ALPHA)
    plt.scatter(y=df["LEfSe-TP"]["WGS-3"], x=df["LEfSe-FP"]["WGS-3"], color="gold", marker='4', alpha=ALPHA)
    plt.scatter(y=df["LEfSe-TP"]["WGS-4"], x=df["LEfSe-FP"]["WGS-4"], color="gold", marker='3', alpha=ALPHA)
    plt.scatter(y=df["LEfSe-TP"]["WGS-5"], x=df["LEfSe-FP"]["WGS-5"], color="gold", marker='2', alpha=ALPHA)

    if MIMIC == True:
        plt.scatter(y=df["TP-miMic"]["Cirrhosis"], x=df["FP-miMic"]["Cirrhosis"], color="hotpink", marker='X',
                    alpha=ALPHA)
        plt.scatter(y=df["TP-miMic"]["Cirrhosis-POP"], x=df["FP-miMic"]["Cirrhosis-POP"], color="hotpink", marker='v',
                    alpha=ALPHA)
        plt.scatter(y=df["TP-miMic"]["ERP020401"], x=df["FP-miMic"]["ERP020401"], color="hotpink", marker='+',
                    alpha=ALPHA)
        plt.scatter(y=df["TP-miMic"]["ERP021216"], x=df["FP-miMic"]["ERP021216"], color="hotpink", marker='D',
                    alpha=ALPHA)
        plt.scatter(y=df["TP-miMic"]["He"], x=df["FP-miMic"]["He"], color="hotpink", marker='d', alpha=ALPHA)
        plt.scatter(y=df["TP-miMic"]["IBD"], x=df["FP-miMic"]["IBD"], color="hotpink", marker='_', alpha=ALPHA)
        plt.scatter(y=df["TP-miMic"]["Jacob"], x=df["FP-miMic"]["Jacob"], color="hotpink", marker='|', alpha=ALPHA)
        plt.scatter(y=df["TP-miMic"]["MF"], x=df["FP-miMic"]["MF"], color="hotpink", marker='x', alpha=ALPHA)
        plt.scatter(y=df["TP-miMic"]["Obesity"], x=df["FP-miMic"]["Obesity"], color="hotpink", marker='^', alpha=ALPHA)
        plt.scatter(y=df["TP-miMic"]["OK30"], x=df["FP-miMic"]["OK30"], color="hotpink", marker='*', alpha=ALPHA)
        plt.scatter(y=df["TP-miMic"]["OK94"], x=df["FP-miMic"]["OK94"], color="hotpink", marker='.', alpha=ALPHA)
        plt.scatter(y=df["TP-miMic"]["PRJNA353587"], x=df["FP-miMic"]["PRJNA353587"], color="hotpink", marker='s',
                    alpha=ALPHA)
        plt.scatter(y=df["TP-miMic"]["PRJNA419097"], x=df["FP-miMic"]["PRJNA419097"], color="hotpink", marker='p',
                    alpha=ALPHA)
        plt.scatter(y=df["TP-miMic"]["T2D"], x=df["FP-miMic"]["T2D"], color="hotpink", marker='8', alpha=ALPHA)
        plt.scatter(y=df["TP-miMic"]["WB"], x=df["FP-miMic"]["WB"], color="hotpink", marker='h', alpha=ALPHA)
        plt.scatter(y=df["TP-miMic"]["WGS-1"], x=df["FP-miMic"]["WGS-1"], color="hotpink", marker='<', alpha=ALPHA)
        plt.scatter(y=df["TP-miMic"]["WGS-2"], x=df["FP-miMic"]["WGS-2"], color="hotpink", marker='>', alpha=ALPHA)
        plt.scatter(y=df["TP-miMic"]["WGS-3"], x=df["FP-miMic"]["WGS-3"], color="hotpink", marker='4', alpha=ALPHA)
        plt.scatter(y=df["TP-miMic"]["WGS-4"], x=df["FP-miMic"]["WGS-4"], color="hotpink", marker='3', alpha=ALPHA)
        plt.scatter(y=df["TP-miMic"]["WGS-5"], x=df["FP-miMic"]["WGS-5"], color="hotpink", marker='2', alpha=ALPHA)

        plt.scatter(y=df["TP-miMic relative"]["Cirrhosis"], x=df["FP-miMic relative"]["Cirrhosis"], color="mediumorchid", marker='X',
                    alpha=ALPHA)
        plt.scatter(y=df["TP-miMic relative"]["Cirrhosis-POP"], x=df["FP-miMic relative"]["Cirrhosis-POP"], color="mediumorchid", marker='v',
                    alpha=ALPHA)
        plt.scatter(y=df["TP-miMic relative"]["ERP020401"], x=df["FP-miMic relative"]["ERP020401"], color="mediumorchid", marker='+',
                    alpha=ALPHA)
        plt.scatter(y=df["TP-miMic relative"]["ERP021216"], x=df["FP-miMic relative"]["ERP021216"], color="mediumorchid", marker='D',
                    alpha=ALPHA)
        plt.scatter(y=df["TP-miMic relative"]["He"], x=df["FP-miMic relative"]["He"], color="mediumorchid", marker='d', alpha=ALPHA)
        plt.scatter(y=df["TP-miMic relative"]["IBD"], x=df["FP-miMic relative"]["IBD"], color="mediumorchid", marker='_', alpha=ALPHA)
        plt.scatter(y=df["TP-miMic relative"]["Jacob"], x=df["FP-miMic relative"]["Jacob"], color="mediumorchid", marker='|', alpha=ALPHA)
        plt.scatter(y=df["TP-miMic relative"]["MF"], x=df["FP-miMic relative"]["MF"], color="mediumorchid", marker='x', alpha=ALPHA)
        plt.scatter(y=df["TP-miMic relative"]["Obesity"], x=df["FP-miMic relative"]["Obesity"], color="mediumorchid", marker='^', alpha=ALPHA)
        plt.scatter(y=df["TP-miMic relative"]["OK30"], x=df["FP-miMic relative"]["OK30"], color="mediumorchid", marker='*', alpha=ALPHA)
        plt.scatter(y=df["TP-miMic relative"]["OK94"], x=df["FP-miMic relative"]["OK94"], color="mediumorchid", marker='.', alpha=ALPHA)
        plt.scatter(y=df["TP-miMic relative"]["PRJNA353587"], x=df["FP-miMic relative"]["PRJNA353587"], color="mediumorchid", marker='s',
                    alpha=ALPHA)
        plt.scatter(y=df["TP-miMic relative"]["PRJNA419097"], x=df["FP-miMic relative"]["PRJNA419097"], color="mediumorchid", marker='p',
                    alpha=ALPHA)
        plt.scatter(y=df["TP-miMic relative"]["T2D"], x=df["FP-miMic relative"]["T2D"], color="mediumorchid", marker='8', alpha=ALPHA)
        plt.scatter(y=df["TP-miMic relative"]["WB"], x=df["FP-miMic relative"]["WB"], color="mediumorchid", marker='h', alpha=ALPHA)
        plt.scatter(y=df["TP-miMic relative"]["WGS-1"], x=df["FP-miMic relative"]["WGS-1"], color="mediumorchid", marker='<', alpha=ALPHA)
        plt.scatter(y=df["TP-miMic relative"]["WGS-2"], x=df["FP-miMic relative"]["WGS-2"], color="mediumorchid", marker='>', alpha=ALPHA)
        plt.scatter(y=df["TP-miMic relative"]["WGS-3"], x=df["FP-miMic relative"]["WGS-3"], color="mediumorchid", marker='4', alpha=ALPHA)
        plt.scatter(y=df["TP-miMic relative"]["WGS-4"], x=df["FP-miMic relative"]["WGS-4"], color="mediumorchid", marker='3', alpha=ALPHA)
        plt.scatter(y=df["TP-miMic relative"]["WGS-5"], x=df["FP-miMic relative"]["WGS-5"], color="mediumorchid", marker='2', alpha=ALPHA)

    plt.plot([df.min().min(), df.max().max()], [df.min().min(), df.max().max()], color='gray', linestyle='--')
    plt.xscale("symlog")
    plt.yscale("symlog")
    plt.xlabel("SP (shuffled positives)",fontsize=SIZE)
    plt.ylabel("RP (real positives)",fontsize=SIZE)
    plt.xticks(fontsize=SIZE)
    plt.yticks(fontsize=SIZE)
    plt.tight_layout()
    if MIMIC:
        plt.savefig("temp/outline_figures/rev1/Fig_1/tp_vs_fp_with_mimic.png")
    else:
        plt.savefig("temp/outline_figures/rev1/Fig_1/tp_vs_fp.png")



    plt.show()
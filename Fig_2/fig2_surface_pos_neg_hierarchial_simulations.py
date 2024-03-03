import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns
from mpl_toolkits.mplot3d import Axes3D
import numpy as np


def plottable_3d_info(df: pd.DataFrame):
    """
    Transform Pandas data into a format that's compatible with
    Matplotlib's surface and wireframe plotting.
    """
    index = df.index
    columns = df.columns

    x, y = np.meshgrid(np.arange(len(columns)), np.arange(len(index)))
    z = np.array([[df[c][i] for c in columns] for i in index])

    xticks = dict(ticks=np.arange(len(columns)), labels=columns)
    yticks = dict(ticks=np.arange(len(index)), labels=index)

    return x, y, z, xticks, yticks


if __name__ == '__main__':
    SIZE = 20
    ST=25
    FS =(5,5)
    TASK = "D"
    mpl.rc('font', family='Times New Roman')
    if TASK == "A":
        fp_1 = pd.read_csv("RESULTS/1D/fp_table_t.csv", index_col=0)
        fp_2 = pd.read_csv("RESULTS/OURS/fp_table_0_0_mu_sis.csv", index_col=0)
        x, y, z, xticks, yticks = plottable_3d_info(fp_1)
        fig = plt.figure(figsize=FS)
        axes = fig.add_subplot(projection='3d')
        s = axes.plot_surface(x, y, z, cmap=plt.get_cmap("Reds"),vmin=0, vmax=25)
        # Adjust xticks to include every third value
        xticks['ticks'] = np.arange(len(xticks['ticks']))[::3]
        xticks['labels'] = xticks['labels'][::3]
        yticks['ticks'] = np.arange(len(yticks['ticks']))[::2]
        yticks['labels'] = yticks['labels'][::2]
        plt.xticks(**xticks,fontsize=SIZE)
        plt.yticks(**yticks,fontsize=SIZE)


        x_f,y_f,z_f,xticks_f, yticks_f= plottable_3d_info(fp_2)
        # Adjust xticks_f to include every third value
        xticks_f['ticks'] = np.arange(len(xticks_f['ticks']))[::3]
        xticks_f['labels'] = xticks_f['labels'][::3]
        s = axes.plot_surface(x_f, y_f, z_f, cmap=plt.get_cmap("Reds"),vmin=0, vmax=25)
        axes.set_xlabel(r'$\mu$',fontsize=SIZE)
        axes.set_ylabel('Number of samples',fontsize=SIZE)
        axes.set_zlabel('FP',fontsize=SIZE)
        # Manually set the coordinates for the text annotations
        leaf_label_coords = (x[0, 0], y[0, 0], z[0, 0]+7)
        mimic_label_coords = (x_f[0, 0], y_f[0, 0], z_f[0, 0]-2)

        # Annotate the surfaces with labels
        axes.text(*leaf_label_coords, "Leaf", fontsize=SIZE, color='black')
        axes.text(*mimic_label_coords, "miMic", fontsize=SIZE, color='black')

        # Set the zticks labels size to SIZE
        for tick in axes.get_zticklabels():
            tick.set_fontsize(SIZE)
        # Add the title
        axes.set_title("Simulation-(0,0,$\mu$) FPs", fontsize=ST,fontweight="bold")
        # Set the azimuth and elevation angles to control the camera view
        azimuth = -13 # Adjust this value to change the azimuth angle (default: 30)
        elevation = 21  # Adjust this value to change the elevation angle (default: 20)


        axes.view_init(elev=elevation, azim=azimuth)
        plt.tight_layout()
        plt.savefig("temp/outline_figures/Fig_2/surface/neg_pos/0_0_mu_fp.png")
        plt.show()


    elif TASK == "B":  # 0,0,mu positives

        tp_2 = pd.read_csv("RESULTS/OURS/tp_table_0_0_mu_sis.csv", index_col=0)

        tp_1 = pd.read_csv("RESULTS/1D/tp_table_t.csv", index_col=0)

        x, y, z, xticks, yticks = plottable_3d_info(tp_1)

        fig = plt.figure(figsize=FS)

        axes = fig.add_subplot(projection='3d')

        s = axes.plot_surface(x, y, z, cmap=plt.get_cmap("Blues"), vmin=0, vmax=100)

        # Adjust xticks to include every third value

        xticks['ticks'] = np.arange(len(xticks['ticks']))[::3]

        xticks['labels'] = xticks['labels'][::3]

        yticks['ticks'] = np.arange(len(yticks['ticks']))[::2]

        yticks['labels'] = yticks['labels'][::2]

        plt.xticks(**xticks, fontsize=SIZE)

        plt.yticks(**yticks, fontsize=SIZE)

        x_f, y_f, z_f, xticks_f, yticks_f = plottable_3d_info(tp_2)

        s = axes.plot_surface(x_f, y_f, z_f, cmap=plt.get_cmap("Blues"), vmin=0, vmax=100)

        axes.set_xlabel(r'$\mu$', fontsize=SIZE)

        axes.set_ylabel('Number of samples', fontsize=SIZE)

        axes.set_zlabel('TP', fontsize=SIZE)

        # Manually set the coordinates for the text annotations

        leaf_label_coords = (x[0, 0], y[0, 0], z[0, 0]+10)

        mimic_label_coords = (x_f[0, 0], y_f[0, 0]-0.5, z_f[0, 0]-15)

        # Annotate the surfaces with labels

        axes.text(*leaf_label_coords, "Leaf", fontsize=SIZE, color='black')

        axes.text(*mimic_label_coords, "miMic", fontsize=SIZE, color='black')

        # Set the zticks labels size to SIZE

        for tick in axes.get_zticklabels():
            tick.set_fontsize(SIZE)

        # Add the title

        axes.set_title(r"Simulation-(0,0,$\mu$) TPs", fontsize=ST, fontweight="bold")

        # Set the azimuth and elevation angles to control the camera view

        azimuth = -13  # Adjust this value to change the azimuth angle (default: 30)

        elevation = 21 # Adjust this value to change the elevation angle (default: 20)

        axes.view_init(elev=elevation, azim=azimuth)

        plt.tight_layout()

        plt.savefig("temp/outline_figures/Fig_2/surface/neg_pos/0_0_mu_tp.png")

        plt.show()

    elif TASK=="C":#
        fp_1 = pd.read_csv("RESULTS/1D/fp_table_relative.csv", index_col=0)
        fp_2 = pd.read_csv("RESULTS/OURS/fp_table_relative_sis.csv", index_col=0)
        x, y, z, xticks, yticks = plottable_3d_info(fp_1)
        fig = plt.figure(figsize=FS)
        axes = fig.add_subplot(projection='3d')
        s = axes.plot_surface(x, y, z, cmap=plt.get_cmap("Reds"), vmin=0, vmax=11)
        # Adjust xticks to include every third value
        xticks['ticks'] = np.arange(len(xticks['ticks']))[::4]
        xticks['labels'] = xticks['labels'][::4]
        yticks['ticks'] = np.arange(len(yticks['ticks']))[::2]
        yticks['labels'] = yticks['labels'][::2]
        plt.xticks(**xticks, fontsize=SIZE)
        plt.yticks(**yticks, fontsize=SIZE)

        x_f, y_f, z_f, xticks_f, yticks_f = plottable_3d_info(fp_2)
        s = axes.plot_surface(x_f, y_f, z_f, cmap=plt.get_cmap("Reds"), vmin=0, vmax=11)
        axes.set_xlabel(r'$\alpha$', fontsize=SIZE)
        axes.set_ylabel('Number of samples', fontsize=SIZE)
        axes.set_zlabel('FP', fontsize=SIZE)
        # Manually set the coordinates for the text annotations

        leaf_label_coords = (x[0, 0], y[0, 0], z[0, 0])

        mimic_label_coords = (x_f[0, 0]+3, y_f[0, 0]-2, z_f[0, 0]+3)

        # Annotate the surfaces with labels

        axes.text(*leaf_label_coords, "Leaf", fontsize=SIZE, color='black')

        axes.text(*mimic_label_coords, "miMic", fontsize=SIZE, color='black')
        # Set the zticks labels size to SIZE
        for tick in axes.get_zticklabels():
            tick.set_fontsize(SIZE)
        # Add the title
        axes.set_title(r"Simulation-(0,$\mu$,$\alpha$*$\mu$) FPs", fontsize=ST, fontweight="bold")
        # Set the azimuth and elevation angles to control the camera view
        azimuth = -14  # Adjust this value to change the azimuth angle (default: 30)
        elevation = 16  # Adjust this value to change the elevation angle (default: 20)
        axes.view_init(elev=elevation, azim=azimuth)
        plt.tight_layout()
        plt.savefig("temp/outline_figures/Fig_2/surface/neg_pos/amu_mu_fp.png")
        plt.show()

    elif TASK=="D":
        tp_2 = pd.read_csv("RESULTS/OURS/tp_table_relative_sis.csv", index_col=0) / 2
        tp_1 = pd.read_csv("RESULTS/1D/tp_table_relative.csv", index_col=0) / 2

        x, y, z, xticks, yticks = plottable_3d_info(tp_1)
        fig = plt.figure(figsize=FS)
        axes = fig.add_subplot(projection='3d')
        s = axes.plot_surface(x, y, z, cmap=plt.get_cmap("Blues"), vmin=0, vmax=100)
        # Adjust xticks to include every third value
        xticks['ticks'] = np.arange(len(xticks['ticks']))[::4]
        xticks['labels'] = xticks['labels'][::4]
        yticks['ticks'] = np.arange(len(yticks['ticks']))[::2]
        yticks['labels'] = yticks['labels'][::2]
        plt.xticks(**xticks, fontsize=SIZE)
        plt.yticks(**yticks, fontsize=SIZE)

        x_f, y_f, z_f, xticks_f, yticks_f = plottable_3d_info(tp_2)
        s = axes.plot_surface(x_f, y_f, z_f, cmap=plt.get_cmap("Blues"),vmin=0, vmax=100)
        axes.set_xlabel(r'$\alpha$', fontsize=SIZE)
        axes.set_ylabel('Number of samples', fontsize=SIZE)
        axes.set_zlabel('TP', fontsize=SIZE)
        # Manually set the coordinates for the text annotations

        leaf_label_coords = (x[0, 0], y[0, 0], z[0, 0]+30)

        mimic_label_coords = (x_f[0, 0], y_f[0, 0], z_f[0, 0])

        # Annotate the surfaces with labels

        axes.text(*leaf_label_coords, "Leaf", fontsize=SIZE, color='black')

        axes.text(*mimic_label_coords, "miMic", fontsize=SIZE, color='black')
        # Set the zticks labels size to SIZE
        for tick in axes.get_zticklabels():
            tick.set_fontsize(SIZE)
        # Add the title
        axes.set_title(r"Simulation-(0,$\mu$,$\alpha$*$\mu$) TPs", fontsize=ST, fontweight="bold")
        azimuth = -14   # Adjust this value to change the azimuth angle (default: 30)
        elevation = 16  # Adjust this value to change the elevation angle (default: 20)
        axes.view_init(elev=elevation, azim=azimuth)
        plt.tight_layout()
        plt.savefig("temp/outline_figures/Fig_2/surface/neg_pos/0_amu_mu_tp.png")
        plt.show()





    c = 0

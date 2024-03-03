import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import statsmodels.api as sm
import numpy as np


if __name__ == '__main__':
    SIZE = 25
    ST =30
    FS =(5,5)
    fs=13
    STATE = "00mu"
    mpl.rc('font', family='Times New Roman')
    if STATE == "000":
        DF0 = pd.read_csv("corr_results/simulations/try_000.csv")

        # Assuming you have a DataFrame df with columns "alpha" and "M"
        X = np.log10(DF0["alpha"])
        Y = np.log10(DF0["M"])

        # Fit the linear regression model
        model = sm.OLS(Y, sm.add_constant(X)).fit()

        # Get the slope (coefficient of "alpha")
        slope = model.params['alpha']

        print("Slope of the regression line:", slope)
        plt.figure(figsize=FS)
        plt.scatter(DF0["alpha"], DF0["M"], color="deeppink")
        plt.xlabel("Leaf confidence", fontsize=SIZE)
        plt.ylabel("miMic confidence", fontsize=SIZE)
        plt.title(r"Regime (0,0,0)", fontsize=ST, fontweight="bold")
        plt.xscale("log")
        plt.yscale("log")
        plt.xticks(fontsize = SIZE)
        plt.yticks(fontsize = SIZE)
        # Plot the regression line
        x_values = np.linspace(min(DF0["alpha"]), max(DF0["alpha"]), 100)
        y_values = slope * np.log10(x_values)
        plt.plot(x_values, 10 ** y_values, color="deeppink", label=f"S = {2.127}")
        plt.legend(fontsize=SIZE)
        plt.tight_layout()
        plt.savefig("corr_results/simulations/new/000.png")
        plt.show()

    elif STATE == "00mu":
        DF1 = pd.read_csv("corr_results/simulations/alpha_005.csv",index_col=0)
        # Assuming you have a DataFrame df with columns "alpha" and "M"
        X = np.log10(DF1["Leaf"])
        Y = np.log10(DF1["miMic"])

        # Fit the linear regression model
        model = sm.OLS(Y, sm.add_constant(X)).fit()

        # Get the slope (coefficient of "alpha")
        slope = model.params["Leaf"]

        print("Slope of the regression line:", slope)
        plt.figure(figsize=FS)
        plt.scatter(DF1["Leaf"], DF1["miMic"], color="deeppink")
        plt.xlabel("Leaf confidence", fontsize=SIZE)
        plt.ylabel("miMic confidence", fontsize=SIZE)
        plt.title(r"Regime (0,0,$\mu$)", fontsize=ST, fontweight="bold")
        plt.xscale("log")
        plt.yscale("log")
        plt.xticks(fontsize=SIZE)
        plt.yticks(fontsize=SIZE)

        # Plot the regression line
        x_values = np.linspace(min(DF1["Leaf"]), max(DF1["Leaf"]), 100)
        y_values = slope * np.log10(x_values)
        plt.plot(x_values, 10 ** y_values, color="deeppink", label=f"S = {1.296}")
        plt.legend(fontsize=SIZE)
        plt.tight_layout()
        plt.savefig("corr_results/simulations/new/00mu.png")
        plt.show()

    elif STATE == "0mumu":

        DF1 = pd.read_csv("corr_results/simulations/new/0mumu_alpha.csv", index_col=0).astype(float)
        # Assuming you have a DataFrame df with columns "alpha" and "M"
        X1 = np.log10(DF1["Leaf"])
        Y1 = np.log10(DF1["miMic-0.25"])
        Y2 = np.log10(DF1["miMic-0.5"])
        Y3 = np.log10(DF1["miMic-1"])
        Y4 = np.log10(DF1["miMic-2"])
        Y5 = np.log10(DF1["miMic-4"])

        # Fit the linear regression model
        model1 = sm.OLS(Y1, sm.add_constant(X1)).fit()
        model2 = sm.OLS(Y2, sm.add_constant(X1)).fit()
        model3 = sm.OLS(Y3, sm.add_constant(X1)).fit()
        model4 = sm.OLS(Y4, sm.add_constant(X1)).fit()
        model5 = sm.OLS(Y5, sm.add_constant(X1)).fit()

        # Get the slope (coefficient of "alpha")
        slope1 = model1.params["Leaf"]
        slope2 = model2.params["Leaf"]
        slope3 = model3.params["Leaf"]
        slope4 = model4.params["Leaf"]
        slope5 = model5.params["Leaf"]

        print("Slope of the regression line:", slope1)
        print("Slope of the regression line:", slope2)
        print("Slope of the regression line:", slope3)
        print("Slope of the regression line:", slope4)
        print("Slope of the regression line:", slope5)
        plt.figure(figsize=FS)
        plt.scatter(DF1["Leaf"], DF1["miMic-0.25"], color="deeppink")
        plt.scatter(DF1["Leaf"], DF1["miMic-0.5"], color="magenta")
        plt.scatter(DF1["Leaf"], DF1["miMic-1"], color="mediumvioletred")
        plt.scatter(DF1["Leaf"], DF1["miMic-2"], color="darkmagenta")
        plt.scatter(DF1["Leaf"], DF1["miMic-4"], color="purple")
        plt.xlabel("Leaf confidence", fontsize=SIZE)
        plt.ylabel("miMic confidence", fontsize=SIZE)
        plt.title(r"Regime (0,$\mu$,$\alpha\mu$)", fontsize=ST, fontweight="bold")
        plt.xscale("log")
        plt.yscale("log")
        plt.xticks(fontsize=SIZE)
        plt.yticks(fontsize=SIZE)

        # Plot the regression line
        x_values = np.linspace(min(DF1["Leaf"]), max(DF1["Leaf"]), 100)
        y_values1 = slope1 * np.log10(x_values)
        y_values2= slope2 * np.log10(x_values)
        y_values3= slope3* np.log10(x_values)
        y_values4= slope4* np.log10(x_values)
        y_values5= slope5* np.log10(x_values)
        plt.plot(x_values, 10 ** y_values1, color="deeppink", label=rf"S(0.25) = {1.342}")
        plt.plot(x_values, 10 ** y_values2, color="magenta", label=rf"S(0.5) = {1.373}")
        plt.plot(x_values, 10 ** y_values3, color="mediumvioletred", label=rf"S(1) = {1.399}")
        plt.plot(x_values, 10 ** y_values4, color="darkmagenta", label=rf"S(2) = {1.405}")
        plt.plot(x_values, 10 ** y_values5, color="purple", label=rf"S(4) = {1.387}")

        plt.legend(fontsize=fs)
        plt.tight_layout()
        plt.savefig("corr_results/simulations/new/0mumu.png")
        plt.show()





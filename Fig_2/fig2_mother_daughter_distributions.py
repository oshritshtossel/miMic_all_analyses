import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

mpl.rc('font', family='Times New Roman')
# Generate data points for the x-axis
x = np.linspace(-5, 5, 1000)
SIZE=25
# Parameters for the two Gaussian distributions
mean1 = 0
std1 = 1

mean2 = 0
std2 = 1 / np.sqrt(3)

# Calculate the probability density function (PDF) for the two Gaussians
pdf1 = (1 / (std1 * np.sqrt(2 * np.pi))) * np.exp(-(x - mean1) ** 2 / (2 * std1 ** 2))
pdf2 = (1 / (std2 * np.sqrt(2 * np.pi))) * np.exp(-(x - mean2) ** 2 / (2 * std2 ** 2))

# Create a plot
plt.figure(figsize=(5, 5))
plt.plot(x, pdf1, label='Daughter\nGaussian',c="deeppink")
plt.plot(x, pdf2, label='Mother\nGaussian',c="mediumvioletred")


# Add labels and legend
plt.xlabel('X',fontsize=SIZE)
plt.ylabel('PDF',fontsize=SIZE)
plt.xticks(fontsize=SIZE)
plt.yticks(fontsize=SIZE)
plt.title('Regime (0,0,0)',fontsize=30,fontweight="bold")
plt.legend(loc="upper left",fontsize=15)#,frameon=False

# Show the plot
plt.grid()
plt.tight_layout()
plt.savefig("temp/outline_figures/Fig_2/gaussians.png")
plt.show()

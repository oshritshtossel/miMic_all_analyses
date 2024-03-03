import matplotlib.pyplot as plt
import matplotlib as mpl
import pandas as pd

SIZE = 15
markers = ['X', 'v', '+', 'D', 'd', '_', '|', 'x', '^', '*', '.', 's', 'p', '8', 'h', '<', '>', '4', '3', '2']

# Set font family to Times New Roman
mpl.rc('font', family='Times New Roman')

# Load data
df = pd.read_csv("Fig_1/Results/all_corrs.csv", index_col=0)
df["significant"] = (df["p-value"] < 0.05)

# Create the bar plot
plt.figure(figsize=(4.8, 4.5))
bars = plt.bar(df.index, df["SCC"], align='center', color="black")

# Manually add asterisks for significant bars
for index, bar in zip(df.index, bars):
    if df["significant"][index]:
        plt.text(bar.get_x() + bar.get_width() / 2, bar.get_height(), '***', ha='center', va='bottom', fontsize=10)

# Plot markers on top of bars
for i, index in enumerate(df.index):
    marker = markers[i % len(markers)]  # Select marker from markers list
    plt.plot(index, df["SCC"][index]/2, marker=marker, markersize=8, color='grey')

# Customize the plot (labels, title, etc.)
plt.ylabel('SCC', fontsize=SIZE)
plt.xticks(rotation=90, fontsize=SIZE)
plt.yticks(fontsize=SIZE)
plt.title(" ", fontsize=SIZE, fontweight="bold")
plt.tight_layout()

# Save the plot
plt.savefig("temp/outline_figures/rev1/Fig_1/sis_corrs.png")
plt.show()

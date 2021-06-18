import pandas as pd 
import numpy as np 
import seaborn as sns
from mpl_toolkits import mplot3d
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.colors import ListedColormap
import matplotlib.pyplot as plt

# Use seaborn themes instead of matplotlib
sns.set()

# ============================================================================ #
filename = '../output/three_biofilm_gates.csv'

# Open the file
df = pd.read_csv(filename)

# Transfer gate labels into numbers
df = df.loc[(df['gate_label'] != "UNKNOWN")]
df.loc[df.gate_label == 'AND', 'gate_label'] = 1
df.loc[df.gate_label == 'OR', 'gate_label'] = 2
df.loc[df.gate_label == 'NOR', 'gate_label'] = 3
df.loc[df.gate_label == 'NAND', 'gate_label'] = 4
df.loc[df.gate_label == 'XOR', 'gate_label'] = 5

cmap = ListedColormap(['green', 'purple', 'red', 'black'])

fig = plt.figure()
ax = fig.add_subplot(111, projection ="3d")
moo = ax.scatter(df['delta'].to_numpy(), 
            df['k'].to_numpy(),
            df['g'].to_numpy(),
            c = df['gate_label'].to_numpy(),
            cmap = cmap)
ax.set_xlabel('Competition Strength')
ax.set_ylabel('Communication Strength')
ax.set_zlabel('Glutamate Concentration')
ax.set_title('Gate Level Description of 3 Biofilm System')

# Show colorbar 
cbar = fig.colorbar(moo, ticks=[2.4, 3.15, 3.9, 4.6])
cbar.ax.set_yticklabels(['OR', 'NOR', 'NAND', 'XOR'])
cbar.ax.set_ylabel('Gate Labels', rotation=270)

plt.show()



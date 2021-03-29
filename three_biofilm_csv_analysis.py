import pandas as pd 
import numpy as np 
import seaborn as sns
from mpl_toolkits import mplot3d
import matplotlib.pyplot as plt

# Use seaborn themes instead of matplotlib
sns.set()

'''
This is the file I will use in order to interpret the data that I gathered in
three_biofilm_analysis.py 
'''

filename = "output/three_biofilm_analysis.csv"

# Read in the file
df = pd.read_csv(filename)

# Separate into 4 dfs depending on input phases
zero_zero = df.loc[(df['biofilm_1_initial_phase'] == 0) & (df['biofilm_3_initial_phase'] == 0)]
zero_one = df.loc[(df['biofilm_1_initial_phase'] == 0) & (df['biofilm_3_initial_phase'] > 3)]
one_zero = df.loc[(df['biofilm_1_initial_phase'] > 3) & (df['biofilm_3_initial_phase'] == 0)]
one_one = df.loc[(df['biofilm_1_initial_phase'] > 3) & (df['biofilm_3_initial_phase'] > 3)]

# Reset the indexes and get rid of extraneous cols
zero_zero = zero_zero.reset_index(drop=True)
zero_one = zero_one.reset_index(drop=True)
one_zero = one_zero.reset_index(drop=True)
one_one = one_one.reset_index(drop=True)

zero_zero = zero_zero.drop(['biofilm_1_initial_phase', 'biofilm_3_initial_phase'], axis=1)
zero_one = zero_one.drop(['biofilm_1_initial_phase', 'biofilm_3_initial_phase'], axis=1)
one_zero = one_zero.drop(['biofilm_1_initial_phase', 'biofilm_3_initial_phase'], axis=1)
one_one = one_one.drop(['biofilm_1_initial_phase', 'biofilm_3_initial_phase'], axis=1)

print(zero_zero)
# print(zero_one)
# print(one_zero)
# print(one_one)

# Zero Zero input (Convert phase differences into boolean indicators)
THRESHOLD_VAL = 1

zero_zero['phase_dif_1_2'] = (zero_zero['phase_dif_1_2'] > THRESHOLD_VAL).astype(int)
zero_zero['phase_dif_1_3'] = (zero_zero['phase_dif_1_3'] > THRESHOLD_VAL).astype(int)
zero_zero['phase_dif_2_3'] = (zero_zero['phase_dif_2_3'] > THRESHOLD_VAL).astype(int)
one_zero['phase_dif_1_2'] = (one_zero['phase_dif_1_2'] > THRESHOLD_VAL).astype(int)
one_zero['phase_dif_1_3'] = (one_zero['phase_dif_1_3'] > THRESHOLD_VAL).astype(int)
one_zero['phase_dif_2_3'] = (one_zero['phase_dif_2_3'] > THRESHOLD_VAL).astype(int)
zero_one['phase_dif_1_2'] = (zero_one['phase_dif_1_2'] > THRESHOLD_VAL).astype(int)
zero_one['phase_dif_1_3'] = (zero_one['phase_dif_1_3'] > THRESHOLD_VAL).astype(int)
zero_one['phase_dif_2_3'] = (zero_one['phase_dif_2_3'] > THRESHOLD_VAL).astype(int)
one_one['phase_dif_1_2'] = (one_one['phase_dif_1_2'] > THRESHOLD_VAL).astype(int)
one_one['phase_dif_1_3'] = (one_one['phase_dif_1_3'] > THRESHOLD_VAL).astype(int)
one_one['phase_dif_2_3'] = (one_one['phase_dif_2_3'] > THRESHOLD_VAL).astype(int)

# ============================================================================ #
fig = plt.figure()

ax = fig.add_subplot(2,2,1, projection ="3d")
ax.scatter(zero_zero['competition_strength'].to_numpy(), 
           zero_zero['communication_strength'].to_numpy(),
           zero_zero['glutamate_concentration'].to_numpy(),
           c = zero_zero['phase_dif_1_2'].to_numpy(),
           cmap = 'binary')
ax.set_xlabel('Competition Strength')
ax.set_ylabel('Communication Strength')
ax.set_zlabel('Glutamate Concentration')
ax.set_title('Zero Zero Input into 3 Biofilm System')


ax = fig.add_subplot(2,2,2, projection ="3d")
ax.scatter(zero_one['competition_strength'].to_numpy(), 
           zero_one['communication_strength'].to_numpy(),
           zero_one['glutamate_concentration'].to_numpy(),
           c = zero_one['phase_dif_1_2'].to_numpy(),
           cmap = 'binary')
ax.set_xlabel('Competition Strength')
ax.set_ylabel('Communication Strength')
ax.set_zlabel('Glutamate Concentration')
ax.set_title('Zero One Input into 3 Biofilm System')


ax = fig.add_subplot(2,2,3, projection ="3d")
ax.scatter(one_zero['competition_strength'].to_numpy(), 
           one_zero['communication_strength'].to_numpy(),
           one_zero['glutamate_concentration'].to_numpy(),
           c = one_zero['phase_dif_1_2'].to_numpy(),
           cmap = 'binary')
ax.set_xlabel('Competition Strength')
ax.set_ylabel('Communication Strength')
ax.set_zlabel('Glutamate Concentration')
ax.set_title('One Zero Input into 3 Biofilm System')


ax = fig.add_subplot(2,2,4, projection ="3d")
ax.scatter(one_one['competition_strength'].to_numpy(), 
           one_one['communication_strength'].to_numpy(),
           one_one['glutamate_concentration'].to_numpy(),
           c = one_one['phase_dif_1_2'].to_numpy(),
           cmap = 'binary')
ax.set_xlabel('Competition Strength')
ax.set_ylabel('Communication Strength')
ax.set_zlabel('Glutamate Concentration')
ax.set_title('One One Input into 3 Biofilm System')

plt.show()




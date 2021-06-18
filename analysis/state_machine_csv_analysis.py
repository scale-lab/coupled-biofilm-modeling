import pandas as pd 
import numpy as np 
import seaborn as sns
from mpl_toolkits import mplot3d
import matplotlib.pyplot as plt

# Use seaborn themes instead of matplotlib
sns.set()

# Read in the file
filename = "../output/three_biofilm_fsm_long_test.csv"
df = pd.read_csv(filename)
print(df)

# g = 0
# k = 0
# delta = 0
# individual_row = df.loc[(df['glutamate_concentration'] == g) & (df['communication_strength'] == k) & (df['competition_strength'] == delta)]
# 

# I think that I should first try and clamp things into one value so I can graph things 
# Let us take the middle as our ground truth. That is, we use it to compare what state we are in
# If phase dif < threshold, in phase; if phase dif > threshold, out of phase
# This means we only have 4 states possible

# if both les
df = df.drop(['mean_output_phase1', 'mean_output_phase2', 'mean_output_phase3', 'phase_dif_1_3'], axis=1)

# Replace biofilm initial phase with easier to read things
df['biofilm_1_initial_phase'] = np.where((df['biofilm_1_initial_phase'] == 0), 0, 1)
df['biofilm_2_initial_phase'] = np.where((df['biofilm_2_initial_phase'] == 0), 0, 1)
df['biofilm_3_initial_phase'] = np.where((df['biofilm_3_initial_phase'] == 0), 0, 1)

THRESHOLD_1_2 = df['phase_dif_1_2'].mean()
THRESHOLD_2_3 = df['phase_dif_2_3'].mean()
mask1 = (df.phase_dif_1_2 > THRESHOLD_1_2) & (df.phase_dif_2_3 < THRESHOLD_2_3)     # 100
mask0 = (df.phase_dif_1_2 < THRESHOLD_1_2) & (df.phase_dif_2_3 < THRESHOLD_2_3)     # 000
mask2 = (df.phase_dif_1_2 < THRESHOLD_1_2) & (df.phase_dif_2_3 > THRESHOLD_2_3)     # 001
mask3 = (df.phase_dif_1_2 > THRESHOLD_1_2) & (df.phase_dif_2_3 > THRESHOLD_2_3)     # 101

df['new'] = np.where(mask0, 0, 'unknown')
df['new'] = np.where(mask1, 1, df.new)
df['new'] = np.where(mask2, 2, df.new)
df['new'] = np.where(mask3, 3, df.new)

# For now, drop all rows that are not 000, 100, 001, or 101
# That is, 111, 110, 011, 010
mask4 = (df.biofilm_1_initial_phase == 1) & (df.biofilm_2_initial_phase == 1) & (df.biofilm_3_initial_phase == 1)
mask5 = (df.biofilm_1_initial_phase == 1) & (df.biofilm_2_initial_phase == 1) & (df.biofilm_3_initial_phase == 0)
mask6 = (df.biofilm_1_initial_phase == 0) & (df.biofilm_2_initial_phase == 1) & (df.biofilm_3_initial_phase == 1)
mask7 = (df.biofilm_1_initial_phase == 0) & (df.biofilm_2_initial_phase == 1) & (df.biofilm_3_initial_phase == 0)
total_mask = mask4 | mask5 | mask6 | mask7

df = df.drop(df[total_mask].index)
print(df)

# Sort values by competition strength and communication strength (express phenotype)
df = df.sort_values(['competition_strength', 'communication_strength'], ascending=True)
print(df)

print('MEAN Difference between Biofilm 1 and 2:', df['phase_dif_1_2'].mean())
print('MEAN Difference between Biofilm 2 and 3:', df['phase_dif_2_3'].mean())

# Output
df.to_csv(r'../output/three_biofilm_fsm_long_test_cleaned_sorted.csv', index=False, header=True)


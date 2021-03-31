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

def find_label(first, second, third, fourth):
    # Compute boolean values
    and_gate = (not first) and (not second) and (not third) and (fourth)
    or_gate = (not first) and second and third and fourth
    nor_gate = first and (not second) and (not third) and (not fourth)
    nand_gate = first and second and third and (not fourth)
    xor_gate = (not first) and second and third and (not fourth)
    unknown = (not and_gate) and (not or_gate) and (not nor_gate) and (not nand_gate) and (not xor_gate)

    # Assign numbers to each one
    if and_gate:
        label = "AND"
    elif or_gate:
        label = "OR"
    elif nor_gate:
        label = "NOR"
    elif nand_gate:
        label = "NAND"
    elif xor_gate:
        label = "XOR"
    else:
        label = "UNKNOWN"

    return label


# ============================================================================ #

filename = "output/three_biofilm_analysis_slower.csv"

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

zero_zero['phase_dif_1_2_bin'] = (zero_zero['phase_dif_1_2'] > THRESHOLD_VAL).astype(int)
zero_zero['phase_dif_1_3_bin'] = (zero_zero['phase_dif_1_3'] > THRESHOLD_VAL).astype(int)
zero_zero['phase_dif_2_3_bin'] = (zero_zero['phase_dif_2_3'] > THRESHOLD_VAL).astype(int)
one_zero['phase_dif_1_2_bin'] = (one_zero['phase_dif_1_2'] > THRESHOLD_VAL).astype(int)
one_zero['phase_dif_1_3_bin'] = (one_zero['phase_dif_1_3'] > THRESHOLD_VAL).astype(int)
one_zero['phase_dif_2_3_bin'] = (one_zero['phase_dif_2_3'] > THRESHOLD_VAL).astype(int)
zero_one['phase_dif_1_2_bin'] = (zero_one['phase_dif_1_2'] > THRESHOLD_VAL).astype(int)
zero_one['phase_dif_1_3_bin'] = (zero_one['phase_dif_1_3'] > THRESHOLD_VAL).astype(int)
zero_one['phase_dif_2_3_bin'] = (zero_one['phase_dif_2_3'] > THRESHOLD_VAL).astype(int)
one_one['phase_dif_1_2_bin'] = (one_one['phase_dif_1_2'] > THRESHOLD_VAL).astype(int)
one_one['phase_dif_1_3_bin'] = (one_one['phase_dif_1_3'] > THRESHOLD_VAL).astype(int)
one_one['phase_dif_2_3_bin'] = (one_one['phase_dif_2_3'] > THRESHOLD_VAL).astype(int)

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

#plt.show()

# ============================================================================ #
# See if I can make it easier to analyze everything, so this makes a new csv
# document that contains information necessary to characterize particular combinations
# as logic gates by.
# FIXME: Right now, we only factor in the phase difference between 1 and 2

OUTPUT_FILE = 'output/three_biofilm_gates.csv'
NUM_ROWS = 250
out = open(OUTPUT_FILE, 'w')
out.write("g,k,delta,o_first,o_second,o_third,o_fourth,first,second,third,fourth,gate_label\n")

# Grab numpy array of g, k, delta
g_vec = zero_zero['glutamate_concentration'].to_numpy()
k_vec = zero_zero['communication_strength'].to_numpy()
delta_vec = zero_zero['competition_strength'].to_numpy()

ofirst_vec = zero_zero['phase_dif_1_2'].to_numpy()
osecond_vec = zero_one['phase_dif_1_2'].to_numpy()
othird_vec = one_zero['phase_dif_1_2'].to_numpy()
ofourth_vec = one_one['phase_dif_1_2'].to_numpy()

first_vec = zero_zero['phase_dif_1_2_bin'].to_numpy()
second_vec = zero_one['phase_dif_1_2_bin'].to_numpy()
third_vec = one_zero['phase_dif_1_2_bin'].to_numpy()
fourth_vec = one_one['phase_dif_1_2_bin'].to_numpy()

# Loop through all rows
for i in range(NUM_ROWS):
    g = g_vec[i]
    k = k_vec[i]
    delta = delta_vec[i]

    ofirst = ofirst_vec[i]
    osecond = osecond_vec[i]
    othird = othird_vec[i]
    ofourth = ofourth_vec[i]

    first = first_vec[i]
    second = second_vec[i]
    third = third_vec[i]
    fourth = fourth_vec[i]

    label = find_label(first, second, third, fourth)

    string = str(g) + "," \
           + str(k) + "," \
           + str(delta) + "," \
           + str(ofirst) + "," \
           + str(osecond) + "," \
           + str(othird) + "," \
           + str(ofourth) + "," \
           + str(first) + "," \
           + str(second) + "," \
           + str(third) + "," \
           + str(fourth) + "," \
           + str(label) + '\n'

    out.write(string)

out.close()




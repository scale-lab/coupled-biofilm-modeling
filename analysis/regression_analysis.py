import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import statsmodels.api as sm
from statsmodels.tools import eval_measures

sns.set(style="ticks")
sns.set(font_scale=1.5)

# ============================================================================ #
def linear_regression(x, y):
    # Try to see if statsmodels has anything different it can provide us
    modified_x = sm.add_constant(x)
    model = sm.OLS(y, modified_x)
    results = model.fit()
    return results

def extract_values(data_df):
    data = pd.DataFrame()
    data["p-value"] = data_df.pvalues
    data["coefficient"] = data_df.params
    print("Phase 1 and 2 difference regression analysis")
    return data


FILENAME = '../output/three_biofilm_fsm_long_test.csv'
df = pd.read_csv(FILENAME)

# Normalize df
df = (df - df.mean()) / df.std()

x = df[['glutamate_concentration','competition_strength','communication_strength','biofilm_1_initial_phase','biofilm_2_initial_phase','biofilm_3_initial_phase']]
# Also might be good to try a non linear model as this may not capture everything
y1 = df['phase_dif_1_2']
y2 = df['phase_dif_2_3']
y3 = df['phase_dif_1_3']

# Perform regression analysis
result1 = linear_regression(x, y1)
result2 = linear_regression(x, y2)
result3 = linear_regression(x, y3)

# Print Regression Analysis to see what our inputs fit
print(result1.summary())
print(result2.summary())
print(result3.summary())

# Extract out the P values and the Coefficients 
values1 = extract_values(result1) 
values2 = extract_values(result2)
values3 = extract_values(result3)

# Print out what we see
print("Phase Dif 1, 2")
print(values1)

print("Phase Dif 2, 3")
print(values2)

print("Phase Dif 1, 3")
print(values3)






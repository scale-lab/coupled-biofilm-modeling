# We will try and extrapolate out a model and see what happens
import numpy as np 
import math 
from scipy.integrate import odeint 
import matplotlib.pyplot as plt
import seaborn as sns 

# Use seaborn themes instead of matplotlib
sns.set()

# Declare global constants (All Basal rates)
alpha = 1.78            # Maximum glutamate consumption rate
w_0 = 1.0               # Basal Intrinsic Frequency
delta_g = 1.2           # Glutamate Consumption Rate
delta_r = 4.46          # Biomass degradation rate
k_r = 1.27              # Saturation threshold for biomass degradation
k_g = 0.01              # Saturation threshold for glutamate consumption
K_0 = 0.016             # Maximum coupling strength FIXME: (This is set at the trkA mutation)
k_theta = 0.29          # Threshold for glutamate modulation of coupling strength
delta_w_0 = 1           # Maximal Glutamate-induced frequency shift
k_omega = 0.19          # Glutamate threshold inhibition of frequency shift
G_t = 1                 # External Glutate Concentration
beta = 6.37                # Glutmate flow rate

def step(x):
    '''
    Definition of step function but for evaluation at a particular value, x
    '''
    if x > 0:
        return x
    else:
        return 0

def kuramato(G):
    '''
    Describes coupling strength for biofilms which is dependent on the amount 
    of current Glutamate in the system
    '''
    os = (K_0 * G) / (k_theta + G)
    return os 

def d_omega(G, theta):
    '''
    Expresses change in the frequency of some arbitrary biofilm
    '''
    cos_out = np.cos(theta)
    denominator = 1 + (G / k_omega)
    out = (delta_w_0 / denominator) * cos_out * step(cos_out)
    return out 

def g_con(G, theta):
    '''
    Refers to glutamate consumption for a particular biofilm, representative of
    theta
    '''
    return ((alpha * G) / (k_g + G)) * (1 - np.sin(theta))

def g_add(g):
    return beta * g

def model(z, t):
    '''
    Model for all of the differential equations where biofilms arranged like:

                [Biofilm 1      Biofilm 2      Biofilm 3]

    Note: z = [theta_1, theta_2, theta_3, G, r1, r2, r3]   
                 0          1       2     3   4   5   6
    '''
    # Compute things I will need multiple times
    kuramato_strength = kuramato(z[3])
    consume_1 = g_con(z[3], z[0])
    consume_2 = g_con(z[3], z[1])
    consume_3 = g_con(z[3], z[2])

    # Phase Change ODEs
    dtheta1_dt = w_0 + d_omega(z[3], z[0]) + (kuramato_strength * np.sin(z[1] - z[0])) #+ (kuramato_strength * np.sin(z[2] - z[0]))
    dtheta2_dt = w_0 + d_omega(z[3], z[1]) + (kuramato_strength * np.sin(z[0] - z[1])) + (kuramato_strength * np.sin(z[2] - z[1]))  # Extrapolating by saying kuramato coupling strength on both sides of this
    dtheta3_dt = w_0 + d_omega(z[3], z[2]) + (kuramato_strength * np.sin(z[1] - z[2])) #+ (kuramato_strength * np.sin(z[0] - z[2]))

    # Change in Concentration
    biofilm_consumption = consume_1 + consume_2 + consume_3
    intrinsic_metabolic_usage = (delta_g * z[4] * z[3]) + (delta_g * z[5] * z[3]) + (delta_g * z[6] * z[3])
    dGdt = g_add(G_t - z[2]) - biofilm_consumption - intrinsic_metabolic_usage

    # Change in Size ODEs
    dr1dt = consume_1 - ((delta_r * z[4]) / (k_r + z[4]))
    dr2dt = consume_2 - ((delta_r * z[5]) / (k_r + z[5]))
    dr3dt = consume_3 - ((delta_r * z[6]) / (k_r + z[6]))
    return [dtheta1_dt, dtheta2_dt, dtheta3_dt, dGdt, dr1dt, dr2dt, dr3dt] 


# ============================================================================ #
'''
For this extrapolation, I am going to try and have three biofilms in a row.
That is, they will all have the same concentration assuming that the microfluidic
device is small enough in which the glutamate concentration will be the same for
all of them

In addition, we will assume that the one in the middle will be influenced by
the coupling of both of the other ones on each side, but the others will NOT
be influenced by the other on the on the other side. That is, biofilm 1 will
not be influenced by biofilm 3 for now.
          [Biofilm 1      Biofilm 2      Biofilm 3]

This model has an issue at the moment, and I am hesitant to extrapolate any further
since we do not really know much. At the moment, we get graphs where the glutamate 
concentration just dives down to 0 pretty quickly, which for some reason makes
our r values negative.

FIXME: Consumption just goes negative for some reason. I believe under further
investigation, this is coming from the glutamate concentration going negative? 
I wonder whether or not I am missing terms in glutamate flow.
'''

# Initial conditions
# [theta_1, theta_2, theta_3, G, r1, r2, r3]
z0 = [0, 0.1, 0.2, G_t, 0, 0, 0]

# Set up time
t = np.linspace(0, 100, num=1000)

# ODE Solve
z = odeint(model, z0, t)

# Graph the phases
# plt.plot(t, z[:,0], label="theta1")
# plt.plot(t, z[:,1], label="theta2")
# plt.plot(t, z[:,2], label="theta3")
mod_phase1 = np.array([x % (2 * np.pi) for x in z[:,0]])
mod_phase2 = np.array([x % (2 * np.pi) for x in z[:,1]])
mod_phase3 = np.array([x % (2 * np.pi) for x in z[:,2]])
plt.plot(t, mod_phase1, label="phase1")
plt.plot(t, mod_phase2, label="phase2")
plt.plot(t, mod_phase3, label="phase3")


# Graph
plt.plot(t, z[:,3], label="G")
plt.plot(t, z[:,4], label="r1")
plt.plot(t, z[:,5], label="r2")
plt.plot(t, z[:,6], label="r3")
plt.legend()
plt.show()

# ============================================================================ #

# Take a look at consumption rates?
con_1 = np.array([g_con(x[0], x[1]) for x in zip(z[:,3], z[:,0])]) 
con_2 = np.array([g_con(x[0], x[1]) for x in zip(z[:,3], z[:,1])]) 
con_3 = np.array([g_con(x[0], x[1]) for x in zip(z[:,3], z[:,2])]) 
total_consumption = con_1 + con_2 + con_3

# Intrinsic Metabolic Rates
metabolic_1 = np.array([delta_g * x[0] * x[1] for x in zip(z[:,4], z[:,3])])
metabolic_2 = np.array([delta_g * x[0] * x[1] for x in zip(z[:,5], z[:,3])])
metabolic_3 = np.array([delta_g * x[0] * x[1] for x in zip(z[:,6], z[:,3])])
total_metabolic = metabolic_1 + metabolic_2 + metabolic_3

# G Addition into Glutamate
g_add = z[:,3] * beta

# Graph all these together
plt.plot(t, total_consumption, label='consumption')
plt.plot(t, total_metabolic, label='total metabolic')
plt.plot(t, total_consumption + total_metabolic, label='combined consumption, metabolism')
plt.plot(t, g_add, label='Addition into model')
plt.legend()
plt.show()

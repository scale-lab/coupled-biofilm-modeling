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
K_0 = 0.073             # Maximum coupling strength 
k_theta = 0.29          # Threshold for glutamate modulation of coupling strength
delta_w_0 = 1           # Maximal Glutamate-induced frequency shift
k_omega = 0.19          # Glutamate threshold inhibition of frequency shift
G_t = 1                 # External Glutate Concentration

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
    cos_out = math.cos(theta)
    denominator = 1 + (G / k_omega)
    out = (delta_w_0 / denominator) * cos_out * step(cos_out)
    return out 

def g_con(G, theta):
    '''
    Refers to glutamate consumption for a particular biofilm, representative of
    theta
    '''
    return ((alpha * G) / (k_g + G)) * (1 - math.sin(theta))

def g_add(g):
    # TODO: Missing equation, emailed authors and waiting on response
    return g

def model(z, t):
    '''
    Model for all of the differential equations. 

    FIXME: I do not know whether or not we actually need a t in there since it
    is not used in any of the equations? Regardless, I will leave it in there
    for now.
    '''
    dtheta1_dt = w_0 + d_omega(z[2], z[0]) + (kuramato(z[2]) * math.sin(z[1] - z[0]))
    dtheta2_dt = w_0 + d_omega(z[2], z[1]) + (kuramato(z[2]) * math.sin(z[0] - z[1]))
    dGdt = g_add(G_t - z[2]) - g_con(z[2], z[0]) - g_con(z[2], z[1]) - (delta_g * z[3] * z[2]) - (delta_g * z[4] * z[2])
    dr1dt = g_con(z[2], z[0]) - ((delta_r * z[3]) / (k_r + z[3]))
    dr2dt = g_con(z[2], z[1]) - ((delta_r * z[4]) / (k_r + z[4]))
    return [dtheta1_dt, dtheta2_dt, dGdt, dr1dt, dr2dt] 


# ============================================================================ #

# Initial conditions
# [theta_1, theta_2, G, r1, r2]
z0 = [0, 0.1, G_t, 0, 0]

# Set up time
t = np.linspace(0, 5)

# ODE Solve
z = odeint(model, z0, t)

# Graph the phases
plt.plot(t, z[:,0])
plt.plot(t, z[:,1])
plt.ylabel("Phase")
plt.xlabel("Time")
plt.title("Phase Change of Coupled Biofilm Model")
plt.show()

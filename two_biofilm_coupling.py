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
G_t = 1.0                 # External Glutate Concentration
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
    Model for all of the differential equations. 
    '''
    dtheta1_dt = w_0 + d_omega(z[2], z[0]) + (kuramato(z[2]) * np.sin(z[1] - z[0]))
    dtheta2_dt = w_0 + d_omega(z[2], z[1]) + (kuramato(z[2]) * np.sin(z[0] - z[1]))
    dGdt = g_add(G_t - z[2]) - g_con(z[2], z[0]) - g_con(z[2], z[1]) - (delta_g * z[3] * z[2]) - (delta_g * z[4] * z[2])
    #dGdt = 0
    dr1dt = g_con(z[2], z[0]) - ((delta_r * z[3]) / (k_r + z[3]))
    dr2dt = g_con(z[2], z[1]) - ((delta_r * z[4]) / (k_r + z[4]))
    return [dtheta1_dt, dtheta2_dt, dGdt, dr1dt, dr2dt] 


# ============================================================================ #
'''
FIXME:

Hmm, not sure if this is actually an issue, but at the moment, with the initial 
conditions below, the phase difference is just approaching 0 all the time. 

This is different from what my intuition should suggest but I can't seem to figure
out what might be going wrong that would bring this problem up. That is, when the
stress of a particular biofilm increases (due to a variety of things, but let
us assume the easiest version in which the only stressor is glutamate concentration

We would assume, according to the paper that in the case that we have low glutamate
concentration they would start to go out of phase in order to share resources, but
it seems that that does not really happen. Instead they both just decrease in size
until glutamate concentration can increase, and then at that moment, they just
keep on going) This makes sense at large concentrations of glutamate, and if we 
keep glutamate concentration constant at some large value, we see that they are in
phase overtime.

Though, with small concentration also with zero change in glutamate, we see that
there is not that much different from if Glutamate initial value is large which is
not what we should expect? I believe that we want them to start to shift out
of phase with low glutamate value, but maybe my time scale is way too large.

I am also looking at the supplemental materials and see that there is a little
bit of information on when they actually grab the information that they do, and it
seems that they grab information at steady state, but it doesn't take long for them
to go into a steady state solution.
'''

# Initial conditions
# [theta_1, theta_2, G, r1, r2]
z0 = [0, 0.1, G_t, 0, 0]

# Set up time
t = np.linspace(0, 300, num=5000)

# ODE Solve
z = odeint(model, z0, t)

# Coupling Strength based off of two oscillators
kuramato = np.array([kuramato(g) for g in z[:,2]])
stress_1 = np.array([1 + np.sin(x) for x in z[:,0]])
stress_2 = np.array([1 + np.sin(x) for x in z[:,1]])
consumption_1 = np.array([g_con(x[0], x[1]) for x in zip(z[:,2], z[:,0])])
consumption_2 = np.array([g_con(x[0], x[1]) for x in zip(z[:,2], z[:,1])])

# Create subplots
fig, ax = plt.subplots(2,2)

# Graph the phases
# plt.plot(t, z[:,0], label="theta1")
# plt.plot(t, z[:,1], label="theta2")
ax[0,0].plot(t, np.abs(z[:,0] - z[:,1]), label="phase difference")
ax[0,0].set_ylabel("Phase Difference")
ax[0,0].set_xlabel("Time")
ax[0,0].set_title("Phase Change of Coupled Biofilm Model")

# Graph coupling strength
ax[0,1].plot(t, kuramato, label="coupling strength")
ax[0,1].set_ylabel("Coupling Strength")
ax[0,1].set_xlabel("Time")
ax[0,1].set_title("Coupling Strength of Two Coupled Biofilm Model")

# Sizes of Coupled Biofilms
ax[1,0].plot(t, z[:,3], label="r1")
ax[1,0].plot(t, z[:,4], label="r2")
ax[1,0].plot(t, z[:,2], label="Glutamate Concentration")
ax[1,0].set_ylabel("Sizes")
ax[1,0].set_xlabel("Time")
ax[1,0].set_title("Sizes of Two Coupled Biofilm Model")
ax[1,0].legend()

# Glutamate consumptions
ax[1,1].plot(t, np.abs(consumption_1 - consumption_2), label="Consumption difference")
# ax[1,1].plot(t, consumption_2, label="Consumption 2")
ax[1,1].set_ylabel("Consumption")
ax[1,1].set_xlabel("Time")
ax[1,1].set_title("Consumption Difference of Two Coupled Biofilm Model")
# ax[1,1].legend()

plt.show()

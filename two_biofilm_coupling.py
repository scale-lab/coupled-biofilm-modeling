from scipy.integrate import odeint 
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns 
import parameters as param
import functions as func

# Use seaborn themes instead of matplotlib
sns.set(font_scale=1.7)

def model(z, t):
    '''
    Model for differential equations for two biofilms coupled together using 
    kuramoto oscillation model with changing coupling constant dependent on the
    glutamate concentration present in the media flowing through the biofilms.

    This is the same model as the Coupling between distant biofilms and emergence
    of nutrient-sharing paper that was released.

    Input: z, a numpy array representing the inputs at a certain time t. It is as
           follows: [theta_1, theta_2, G, r1, r2] where theta represents the 
           phase variables for their respective biofilms, G represents the 
           glutamate concentration, and r1, r2 represent the size of the biofilms.
    Input: t, a float representing the time.
    '''
    # Phase Variables
    dtheta1_dt = param.w_0 + func.d_omega(z[2], z[0]) + (func.kuramoto(z[2], None) * np.sin(z[1] - z[0]))
    dtheta2_dt = param.w_0 + func.d_omega(z[2], z[1]) + (func.kuramoto(z[2], None) * np.sin(z[0] - z[1]))

    # Glutamate Consumption Rates
    biomass_consump_1 = func.g_con(z[2], z[0])
    biomass_consump_2 = func.g_con(z[2], z[1])
    total_biomass_consump = biomass_consump_1 + biomass_consump_2
    total_metabolic_consump = func.m_consump(z[3], z[2], None) + func.m_consump(z[4], z[2], None)
    dGdt = func.g_add(param.G_t - z[2]) - total_biomass_consump - total_metabolic_consump

    # Biofilm Growth Rates
    dr1dt = biomass_consump_1 - ((param.delta_r * z[3]) / (param.k_r + z[3]))
    dr2dt = biomass_consump_2 - ((param.delta_r * z[4]) / (param.k_r + z[4]))
    return [dtheta1_dt, dtheta2_dt, dGdt, dr1dt, dr2dt] 

# ============================================================================ #

# Initial conditions
# [theta_1, theta_2, G, r1, r2]
z0 = [0, 0.1, param.G_t, 0, 0]

# Set up time
t = np.linspace(0, 300, num=5000)

# ODE Solve
z = odeint(model, z0, t)

# Creating vectors for various things we may want to graph
kuramoto = np.array([func.kuramoto(g, None) for g in z[:,2]])
# stress_1 = np.array([1 + np.sin(x) for x in z[:,0]])
# stress_2 = np.array([1 + np.sin(x) for x in z[:,1]])
consumption_1 = np.array([func.g_con(x[0], x[1]) for x in zip(z[:,2], z[:,0])])
consumption_2 = np.array([func.g_con(x[0], x[1]) for x in zip(z[:,2], z[:,1])])

# Create subplots
fig, ax = plt.subplots(2,2)

# Graph the phases
# plt.plot(t, z[:,0], label="theta1")
# plt.plot(t, z[:,1], label="theta2")
ax[0,0].plot(t, np.abs(z[:,0] - z[:,1]), label="phase difference")
ax[0,0].set_ylabel("Phase Difference (rad)")
ax[0,0].set_xlabel("Time")
ax[0,0].set_title("Phase Change of Coupled Biofilm Model")

# Graph coupling strength
ax[0,1].plot(t, kuramoto, label="coupling strength")
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
ax[1,1].set_ylabel("Consumption")
ax[1,1].set_xlabel("Time")
ax[1,1].set_title("Consumption Difference of Two Coupled Biofilm Model")

plt.show()

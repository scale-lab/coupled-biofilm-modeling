from scipy.integrate import odeint 
import matplotlib.pyplot as plt
import seaborn as sns 
import numpy as np
import parameters as param
import functions as func

# Use seaborn themes instead of matplotlib
sns.set(font_scale=1.7)

def model(z, t):
    '''
    Model for all of the differential equations where biofilms arranged like:

                [Biofilm 1      Biofilm 2      Biofilm 3]

    Note: z = [theta_1, theta_2, theta_3, G, r1, r2, r3]   
                 0          1       2     3   4   5   6
    '''
    # Compute things I will need multiple times
    kuramoto_strength = func.kuramoto(z[3], None)
    consume_1 = func.g_con(z[3], z[0])
    consume_2 = func.g_con(z[3], z[1])
    consume_3 = func.g_con(z[3], z[2])
    metabolic_1 = func.m_consump(z[4], z[3], None)
    metabolic_2 = func.m_consump(z[5], z[3], None)
    metabolic_3 = func.m_consump(z[6], z[3], None)

    # Phase Change ODEs
    dtheta1_dt = param.w_0 + func.d_omega(z[3], z[0]) + (kuramoto_strength * np.sin(z[1] - z[0])) 
    dtheta2_dt = param.w_0 + func.d_omega(z[3], z[1]) + (kuramoto_strength * np.sin(z[0] - z[1])) + (kuramoto_strength * np.sin(z[2] - z[1]))  # Extrapolating by saying kuramoto coupling strength on both sides of this
    dtheta3_dt = param.w_0 + func.d_omega(z[3], z[2]) + (kuramoto_strength * np.sin(z[1] - z[2])) 

    # Change in Concentration
    biofilm_consumption = consume_1 + consume_2 + consume_3
    intrinsic_metabolic_usage = metabolic_1 + metabolic_2 + metabolic_3
    dGdt = func.g_add(param.G_t - z[3]) - biofilm_consumption - intrinsic_metabolic_usage

    # Change in Size ODEs
    dr1dt = consume_1 - ((param.delta_r * z[4]) / (param.k_r + z[4]))
    dr2dt = consume_2 - ((param.delta_r * z[5]) / (param.k_r + z[5]))
    dr3dt = consume_3 - ((param.delta_r * z[6]) / (param.k_r + z[6]))
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
theta1 = np.pi
theta2 = np.pi
theta3 = np.pi
glutamate = param.G_t
r1 = 0
r2 = 0
r3 = 0
z0 = [theta1, theta2, theta3, glutamate, r1, r2, r3]

# Set up time
t = np.linspace(0, 2000, num=4000)

# ODE Solve
z = odeint(model, z0, t)

# Graph the phase difference between z1, z2 and z1, z3
modulo_phase_1 = np.mod(z[:,0], 2 * np.pi)
modulo_phase_2 = np.mod(z[:,1], 2 * np.pi)
modulo_phase_3 = np.mod(z[:,2], 2 * np.pi)

plt.plot(t, modulo_phase_1, label="mod1")
plt.plot(t, modulo_phase_2, label="mod2")
plt.plot(t, modulo_phase_3, label="mod3")
plt.xlabel("time")
plt.ylabel("Phase Modulo")
plt.title("Phase difference between 3 biofilms over time")
plt.legend()
plt.show()



phase_dif_1_3 = np.abs(modulo_phase_1 - modulo_phase_3)         
phase_dif_1_2 = np.abs(modulo_phase_1 - modulo_phase_2)
phase_dif_2_3 = np.abs(modulo_phase_2 - modulo_phase_3)

plt.plot(t, phase_dif_1_2, label="abs(phase1 - phase2)")
plt.plot(t, phase_dif_2_3, label="abs(phase3 - phase2)")
plt.plot(t, phase_dif_1_3, label="abs(phase1 - phase3)")
plt.xlabel("time")
plt.ylabel("Phase difference (rad)")
plt.title("Phase difference between 3 biofilms over time")
plt.legend()

# Quick and dirty math to try and figure out end of phase
print("Phase 1 Mean:", np.mean(modulo_phase_1[-100:]))
print("Phase 2 Mean:", np.mean(modulo_phase_2[-100:]))
print("Phase 3 Mean:", np.mean(modulo_phase_3[-100:]))
print("phase_dif 1 3 mean:", np.mean(phase_dif_1_3[-100:]))
print("phase_dif 2 3 mean:", np.mean(phase_dif_2_3[-100:]))
print("phase_dif 1 2 mean:", np.mean(phase_dif_1_2[-100:]))



plt.show()

# ============================================================================ #

# Graph of regular variables
# plt.plot(t, z[:,3], label="G")
# plt.plot(t, z[:,4], label="r1")
# plt.plot(t, z[:,5], label="r2")
# plt.plot(t, z[:,6], label="r3")
# plt.legend()
# plt.show()

# ============================================================================ #

# # Take a look at consumption rates?
# con_1 = np.array([func.g_con(x[0], x[1]) for x in zip(z[:,3], z[:,0])]) 
# con_2 = np.array([func.g_con(x[0], x[1]) for x in zip(z[:,3], z[:,1])]) 
# con_3 = np.array([func.g_con(x[0], x[1]) for x in zip(z[:,3], z[:,2])]) 
# total_consumption = con_1 + con_2 + con_3

# # Intrinsic Metabolic Rates
# metabolic_1 = np.array([param.delta_g * x[0] * x[1] for x in zip(z[:,4], z[:,3])])
# metabolic_2 = np.array([param.delta_g * x[0] * x[1] for x in zip(z[:,5], z[:,3])])
# metabolic_3 = np.array([param.delta_g * x[0] * x[1] for x in zip(z[:,6], z[:,3])])
# total_metabolic = metabolic_1 + metabolic_2 + metabolic_3

# # G Addition into Glutamate
# g_add = z[:,3] * param.beta

# # Graph all these together
# plt.plot(t, total_consumption, label='consumption')
# plt.plot(t, total_metabolic, label='total metabolic')
# plt.plot(t, total_consumption + total_metabolic, label='combined consumption, metabolism')
# plt.plot(t, g_add, label='Addition into model')
# plt.legend()
# plt.show()

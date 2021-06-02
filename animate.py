import matplotlib.animation as animation
from scipy.integrate import odeint 
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
from numpy import sin, cos

# Custom Python Functions
import parameters as param
import functions as func

# Animation Parameters
SIZE = 3000
C = cm.twilight_shifted

# Model Definition
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

# Initial conditions
theta1 = np.pi
theta2 = np.pi
theta3 = 0
glutamate = param.G_t
r1 = 0
r2 = 0
r3 = 0
z0 = [theta1, theta2, theta3, glutamate, r1, r2, r3]
num_points = 100
end_time = 50
dt = end_time / num_points
t = np.linspace(0, end_time, num=num_points)

# ODE Solve
z = odeint(model, z0, t)

# Break out output 2D array back into 1D vectors
t1 = z[:, 0]
t2 = z[:, 1]
t3 = z[:, 2]
g = z[:, 3]
size1 = z[:, 4]
size2 = z[:, 5]
size3 = z[:, 6]

# Convert the stress phases into 
modulo_phase_1 = np.mod(z[:,0], 2 * np.pi)
modulo_phase_2 = np.mod(z[:,1], 2 * np.pi)
modulo_phase_3 = np.mod(z[:,2], 2 * np.pi)

# Configure the plot
fig = plt.figure()

# Text templates and initial declarations (where they exist on the map)
time_template = 'time = %.1f'
theta1_template = 'Theta1 =%.2f'
theta2_template = 'Theta2 =%.2f'
theta3_template = 'Theta3 =%.2f'

#===========
def animate(i):
    fig.clear()
    ax = fig.add_subplot(111, autoscale_on=False, xlim=(-2, 2), ylim=(-2, 2))
    
    # Create Titles axes and things
    ax.set_title("MOO")
    ax.set_xticks([])
    ax.set_yticks([])
    
    # Create 1 value scatter plots representing each of the points
    plt.scatter (-1, 0, c = modulo_phase_1[i], vmin = 0, vmax = 2*np.pi, cmap = C, s=SIZE)
    plt.scatter (0, 0, c = modulo_phase_2[i], vmin = 0, vmax = 2*np.pi, cmap = C, s=SIZE)
    plt.scatter (1, 0, c = modulo_phase_3[i], vmin = 0, vmax = 2*np.pi, cmap = C, s=SIZE)

    # Update text in the animation
    ax.text(-1.8, 1.75, time_template % (i*dt))
    ax.text(-1.8, 1.5, theta1_template % modulo_phase_1[i])
    ax.text(-1.8, 1.25, theta2_template % modulo_phase_2[i])
    ax.text(-1.8, 1.0, theta3_template % modulo_phase_3[i])

ani = animation.FuncAnimation(fig, animate, len(t), interval=5)

ani.save('output/three_biofilm.gif', writer=animation.PillowWriter(fps=15))

plt.show()

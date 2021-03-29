from scipy.integrate import odeint 
import numpy as np
import parameters as param
import functions as func
import os

# FIXME: it might be good to change the name of this later so it is easier to 
#        tell what this file is.

'''
This is just a quick file that takes the information from three_biofilm_coupling
and then applies it to a variety of different experiments. That is, we change
the competition strength, coupling strength, and glutamate concentration of 
the external solution and then put this all into a .csv file for further
analysis that will be put in the output folder. 

At the moment, I am just going to vary biofilm 1 and biofilm 3 and say keep the
biofilm 2 phase constant initially at 0 out of phase. While this may not be the
case always, for now I will assume that we will be able to change phases arbitrarily 
which may not be true, but will use this to further our understanding of this 
model
'''

def model(z, t, c_g, c_delta, c_k):
    '''
    c_g, c_delta, c_k refer to the custom glutamate, custom delta_g, custom k_0
    that we define and we change.

    Model for all of the differential equations where biofilms arranged like:

                [Biofilm 1      Biofilm 2      Biofilm 3]

    Note: z = [theta_1, theta_2, theta_3, G, r1, r2, r3]   
                 0          1       2     3   4   5   6

    Currently, this uses nearest neighbors for the coupling model, but we can
    later change that if it really matters
    '''
    # Compute things I will need multiple times
    kuramoto_strength = func.kuramoto(z[3], c_k)
    consume_1 = func.g_con(z[3], z[0])
    consume_2 = func.g_con(z[3], z[1])
    consume_3 = func.g_con(z[3], z[2])
    metabolic_1 = func.m_consump(z[4], z[3], c_delta)
    metabolic_2 = func.m_consump(z[5], z[3], c_delta)
    metabolic_3 = func.m_consump(z[6], z[3], c_delta)

    # Phase Change ODEs
    dtheta1_dt = param.w_0 + func.d_omega(z[3], z[0]) + (kuramoto_strength * np.sin(z[1] - z[0])) 
    dtheta2_dt = param.w_0 + func.d_omega(z[3], z[1]) + (kuramoto_strength * np.sin(z[0] - z[1])) + (kuramoto_strength * np.sin(z[2] - z[1])) 
    dtheta3_dt = param.w_0 + func.d_omega(z[3], z[2]) + (kuramoto_strength * np.sin(z[1] - z[2])) 

    # Change in Concentration
    biofilm_consumption = consume_1 + consume_2 + consume_3
    intrinsic_metabolic_usage = metabolic_1 + metabolic_2 + metabolic_3
    dGdt = func.g_add(c_g - z[3]) - biofilm_consumption - intrinsic_metabolic_usage

    # Change in Size ODEs
    dr1dt = consume_1 - ((param.delta_r * z[4]) / (param.k_r + z[4]))
    dr2dt = consume_2 - ((param.delta_r * z[5]) / (param.k_r + z[5]))
    dr3dt = consume_3 - ((param.delta_r * z[6]) / (param.k_r + z[6]))
    return [dtheta1_dt, dtheta2_dt, dtheta3_dt, dGdt, dr1dt, dr2dt, dr3dt] 

# ============================================================================ #
# Hyperparameters for our testing
TIME_VECTOR_LENGTH = 1000
FINAL_TIME = 1000
AMOUNT_MEAN_POINTS = 10
R1 = 0
R2 = 0
R3 = 0
THETA_2 = 0
G_RANGE = (0.5, 2.5)
K_RANGE = (0.2, 1.2) # Note: these are normalized values based on WT
DELTA_RANGE = (1, 100) # Note: these are normalized values based on WT at 1
G_STEPS = 10
K_STEPS = 5
DELTA_STEPS = 5
OUTPUT_FILE = "output/three_biofilm_analysis_slower.csv" 
TABLE_LABELS = "glutamate_concentration,competition_strength,communication_strength,biofilm_1_initial_phase,biofilm_3_initial_phase,mean_output_phase1,mean_output_phase2,mean_output_phase3,phase_dif_1_2,phase_dif_2_3,phase_dif_1_3\n"

# ============================================================================ #
# Open file descriptor and write first line
f = open(OUTPUT_FILE, "w")
f.write(TABLE_LABELS)

# First calculate the vectors in which we will be using
g_vector = np.linspace(G_RANGE[0], G_RANGE[1], num=G_STEPS)
k_vector = np.linspace(K_RANGE[0], K_RANGE[1], num=K_STEPS)
delta_vector = np.linspace(DELTA_RANGE[0], DELTA_RANGE[1], num=DELTA_STEPS)
phase_vector = [(0, 0), (0, np.pi), (np.pi, 0), (np.pi, np.pi)]

# Setup time vector
t = np.linspace(0, FINAL_TIME, num=TIME_VECTOR_LENGTH)

# Start!
print("Start going through")
# Start looping through all possible combinations
# FIXME: I know it is possible to get rid of this O(n^4) loop using numpy and
#        matrix manipulation to do everything, but this shouldn't take that long
#        to run anyways
for g in g_vector:
    for k in k_vector:
        for delta in delta_vector:
            for p1, p3 in phase_vector:
                # Get rid of normalizations
                real_k = k * param.K_0
                real_delta = delta * param.delta_g

                # Initial conditions
                z0 =  [p1, THETA_2, p3, g, R1, R2, R3]

                # ODE Solve
                z = odeint(model, z0, t, args=(g, real_delta, real_k))

                # Gather the points that we need
                modulo_phase_1 = np.mod(z[(-1 * AMOUNT_MEAN_POINTS):,0], 2 * np.pi)
                modulo_phase_2 = np.mod(z[(-1 * AMOUNT_MEAN_POINTS):,1], 2 * np.pi)
                modulo_phase_3 = np.mod(z[(-1 * AMOUNT_MEAN_POINTS):,2], 2 * np.pi)

                mean_phase_1 = np.mean(modulo_phase_1)
                mean_phase_2 = np.mean(modulo_phase_2)
                mean_phase_3 = np.mean(modulo_phase_3)

                phase_dif_1_2 = np.mean(np.abs(modulo_phase_1 - modulo_phase_2))
                phase_dif_2_3 = np.mean(np.abs(modulo_phase_2 - modulo_phase_3))
                phase_dif_1_3 = np.mean(np.abs(modulo_phase_1 - modulo_phase_3))

                # Construct .csv file input
                line = str(g) + "," \
                     + str(delta) + "," \
                     + str(k) + "," \
                     + str(p1) + "," \
                     + str(p3) + "," \
                     + str(mean_phase_1) + "," \
                     + str(mean_phase_2) + "," \
                     + str(mean_phase_3) + "," \
                     + str(phase_dif_1_2) + "," \
                     + str(phase_dif_2_3) + "," \
                     + str(phase_dif_1_3) + "\n" \
                
                # Write into file
                f.write(line)

            # Display where we are in the console
            print("G: {}, K: {}, Delta: {}".format(g, k, delta))

        # Force file sync so that I don't lose a bunch of information if something happens
        f.flush()
        os.fsync(f.fileno())

# Close file descriptor
print("Finished!")
f.close()
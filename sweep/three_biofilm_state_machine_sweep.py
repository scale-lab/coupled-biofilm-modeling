from scipy.integrate import odeint 
import threading
import numpy as np
import parameters as param
import functions as func
import os

'''
This is just a quick file that takes the information from three_biofilm_coupling
and then applies it to a variety of different experiments. Particularly, this is
a file that is very similar to the three_biofilm_testing.py file, but instead
we try to model biofilms as a finite state machine in which we have some initial
condition and then we go through all of the possible combinations at all different
types of g, delta, and k to see whether or not we can make any useful state 
diagrams out of our biofilm system.

Note: One primary difference is that I am going to try and change the script so
that one can run this on multiple threads so that it takes a lot shorter of 
time to do anything. The number of threads we spawn will kind of be bad code,
but it will spawn based on the number of steps in the glutamate that we have.

Note: This multithreading is at the consequence of having the CSV not be in order.
That is, entries will not be in order, which is not really an issue since we 
will not be looking at the CSV very closely anyways.

Note: According to the ODEINT documentation, it is not possible to do multithreading
with ODEINT as it only instantiates one ODE solver to use (which is actually super
annoying, but I presume I will just use the lab servers when those come back up)
'''

# ============================================================================ #
# Hyperparameters for our testing
TIME_VECTOR_LENGTH = 4000
FINAL_TIME = 2000
AMOUNT_MEAN_POINTS = 100
R1 = 0
R2 = 0
R3 = 0
G_RANGE = (0.5, 2.5)
K_RANGE = (0.2, 1.2) # Note: these are normalized values based on WT
DELTA_RANGE = (1, 100) # Note: these are normalized values based on WT at 1
G_STEPS = 15
K_STEPS = 15
DELTA_STEPS = 15
OUTPUT_FILE = "output/three_biofilm_fsm_long_test.csv" 
TABLE_LABELS = "glutamate_concentration,competition_strength,communication_strength,biofilm_1_initial_phase,biofilm_2_initial_phase,biofilm_3_initial_phase,mean_output_phase1,mean_output_phase2,mean_output_phase3,phase_dif_1_2,phase_dif_2_3,phase_dif_1_3\n"

# First calculate the vectors in which we will be using
TIME_VECTOR = np.linspace(0, FINAL_TIME, num=TIME_VECTOR_LENGTH)
G_VECTOR = np.linspace(G_RANGE[0], G_RANGE[1], num=G_STEPS)
K_VECTOR = np.linspace(K_RANGE[0], K_RANGE[1], num=K_STEPS)
DELTA_VECTOR = np.linspace(DELTA_RANGE[0], DELTA_RANGE[1], num=DELTA_STEPS)
PHASE_VECTOR = [(0, 0, 0), (0, 0, np.pi), (0, np.pi, 0), (0, np.pi, np.pi), (np.pi, 0, 0), (np.pi, 0, np.pi), (np.pi, np.pi, 0), (np.pi, np.pi, np.pi)]

# # Global Lock for writing into list and list with lines to write out
# mutex = threading.Lock()
# lines = []

# ============================================================================ #

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
# Various functions

def compute_odes(f):
    '''
    This is the function that we will pass to each of our threads when we spawn
    them in our main process
    '''
    for g in G_VECTOR:
        for k in K_VECTOR:
            for delta in DELTA_VECTOR:
                for p1, p2, p3 in PHASE_VECTOR:
                    # Get rid of normalizations
                    real_k = k * param.K_0
                    real_delta = delta * param.delta_g

                    # Initial conditions
                    z0 =  [p1, p2, p3, g, R1, R2, R3]

                    # ODE Solve
                    z = odeint(model, z0, TIME_VECTOR, args=(g, real_delta, real_k))

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
                        + str(p2) + "," \
                        + str(p3) + "," \
                        + str(mean_phase_1) + "," \
                        + str(mean_phase_2) + "," \
                        + str(mean_phase_3) + "," \
                        + str(phase_dif_1_2) + "," \
                        + str(phase_dif_2_3) + "," \
                        + str(phase_dif_1_3) + "\n" \
                    
                    # # Write into file
                    # lines.append(line)

                    # Write into file
                    f.write(line)

                # Display where we are in the console
                print("G: {}, K: {}, Delta: {}".format(g, k, delta))

            # Force file sync so that I don't lose a bunch of information if something happens
            f.flush()
            os.fsync(f.fileno())

    # Return after we finish everything
    return

# ============================================================================ #
# Main process
with open(OUTPUT_FILE, "w") as f:
    f.write(TABLE_LABELS)
    compute_odes(f)

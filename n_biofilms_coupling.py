from scipy.integrate import odeint 
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns 
import parameters as param
import functions as func
from k_coloring import txttoMatrix

'''
Some cool notes. I am not sure why I did not just put all of these new
function inside of the functions folder... That would make the most sense
but here we are.

TODO: Move these over to the functions folder.
'''

'''
parse_adjacency_matrix is a function that takes in an numpy array-like
and returns a hash map of lists where the keys refer to the row that
some biofilm is a part of and the list refers to the index of the biofilm
in the graph connection.

We can also handle weighted adjacency matrices, meaning that anything in
the matrix corresponds to the weight, as long as it is not zero.

Input: arr, a weighted adjacency matrix, representing the graph connections
       and couplings
Output: relations, a dictionary of lists, where the key represents the particular
        biofilm and the list represents a bunch of tuples representing 
        (index biofilm connection, weight).
'''
def parse_adjacency_matrix(arr):
    relations = {}
    arr_length = arr.shape[0]

    # Loop through Rows
    for i in range(arr_length):
        i_relations = []

        # For each row, look at columns to gather couplings
        for j in range(arr_length):
            if (arr[i][j] != 0):
                weight = arr[i][j]
                i_relations.append((j, weight))
            
        relations[i] = i_relations

    return relations

'''
kuramoto_coupling is a function that takes in the k_strength, relationship 
dictionary, a vector of state variables in the differential equations, and
the number of biofilms. The purpose of this function is to generate a list of
all of the equations responsible for the differential equation controlling the
phase of a particular biofilm. 

This means that if we have N biofilms, we will have an output that has N equations
corresponding to it. 

Input: k_strength, a float representing the coupling strength based off of 
       hyperparameters. 
Input: relations, a dictionary of lists where each element of the list is a tuple
       represented by (index of biofilm connection, weight). This is often obtained
       by a call to parse_adjacency_matrix()
Input: z, an (2 * num_biofilms) + 1 length vector which includes num_biofilms amount of
       phase variables, num_biofilms amounts of size variables, and one global variable
       on glutamate concentration in the system as a whole.
Input: num_biofilms, an integer that represents the number of biofilms in the whole system
       that we are dealing with.

Output: a list of coupling functions that is num_biofilms length long. Each of the elements
        of the list describe the phase relations between a particular biofilm, i, and its
        connected biofilms.

NOTE:  This is an extrapolation of original biofilm model where we just add all couplings 
       together regardless of distance or anything
'''
def kuramoto_coupling(k_strength, relations, z, num_biofilms):
    couplings = []
    g_index = (num_biofilms*2)

    for i in range(num_biofilms):
        # Get constants first
        dtheta_i = param.w_0 + func.d_omega(z[g_index], z[i])


        for j, weight in relations[i]:
            # i and j coupling
            couple_strength = k_strength * weight
            dtheta_i += (couple_strength * np.sin(z[j] - z[i]))
        
        couplings.append(dtheta_i)

    return couplings

'''
concentration_change is a function that creates the differential equation
responsible for the concentration change based off of glutamate consumption
by cells (metabolism), flow, and consumption due to increasing cell size.

Input: z, state array discussed above. Needed for scipy.Integrate
Input: num_biofilms, the number of biofilms in a particular system.
Output: a list of one, representing the dGdt function. We note that it is
        in a list to make all the differential equation function generators
        consistent so that we can just concatenate everything later.
'''
def concentration_change(z, num_biofilms):
    g_index = (num_biofilms * 2)
    start_of_r = num_biofilms
    dGdt = func.g_add(param.G_t - z[g_index])

    # Loop through and establish metabolic and consumption of glutamate
    for i in range(num_biofilms):
        r_index = start_of_r + i

        consume_i = func.g_con(z[g_index], z[i])
        metabolic_i = func.m_consump(z[r_index], z[g_index], None)

        dGdt -= (consume_i - metabolic_i)
    
    return [dGdt]

'''
size_change is the function that generates all of the differential
equations related to the change in radius of different biofilms.

Input: z, the state vector we had discussed beforehand
Input: num_biofilms, the number of biofilms in the particular system
Output: drdt_list, a list of differential equations representing the
        radius equations.
'''
def size_change(z, num_biofilms):
    start_of_r = num_biofilms
    g_index = (num_biofilms * 2)
    drdt_list = []
    
    # Loop through each biofilm and establish the equation
    for i in range(num_biofilms):
        r_index = start_of_r + i
        consume_i = func.g_con(z[g_index], z[i])

        drdt_i = consume_i - (param.delta_r * z[r_index]) / (param.k_r + z[r_index])

        drdt_list.append(drdt_i)
    
    return drdt_list


'''
create_model is our bread and butter function here. It takes in an adjacency
matrix, parses it, and then establishes the number of differential equations
needed for this system and creates the model that has the attributes of the
weighted adjacency matrix. 

Input: matrix, an np array-like representing a square weighted adjacency matrix
Output: model, a function that describes the differential equations set by the
        adjacency matrix in the context of our biofilm dynamical model.

FIXME: At the moment, our model can only support one type of bacteria at a
       time, meaning that we cannot have different coupling strengths and
       such at the moment. It would not be that hard to fix, but that is 
       just something to consider.
'''
# Note: Z vector as follows: (n biofilm phases, n biofilm radius, glutamate_concentration)
def create_model(matrix):
    # Enforce Square Matrix
    if (matrix.shape[0] != matrix.shape[1]):
        raise Exception("Matrix not Square")

    # Take numpy matrix and produce relations
    relations = parse_adjacency_matrix(matrix)
    num_biofilms = matrix.shape[0]

    # Establish Size of Z that I need (n biofilm phases, n biofilm , gl)
    z_size = 2 * num_biofilms + 1

    # Create Model 
    def model(z,t):
        # TODO: Can replace None with try case for custom params that can be passed in from main func
        couple_strength = func.kuramoto(z[z_size-1], None)

        # Create Phase Change ODE:
        phase_change_odes = kuramoto_coupling(couple_strength, relations, z, num_biofilms)

        # Change in Size ODES
        r_change_odes = size_change(z, num_biofilms)

        # Create concentration change ODE
        c_change_ode = concentration_change(z, num_biofilms)

        all_odes = phase_change_odes + r_change_odes + c_change_ode
        
        return all_odes

    return model
        

# =======================================

# Some Test Code for the functions

z0 = []
t = np.linspace(0, 400, num=2000)

filename = "matrix.txt"
matrix = np.asarray(txttoMatrix(filename))
num_biofilms = matrix.shape[0]

# Sample Initial Conditions

for i in range(num_biofilms):
    theta_i = 0
    z0.append(theta_i)

z0[0] = 0
z0[1] = np.pi
z0[2] = 0

for j in range(num_biofilms):
    rad_i = 0
    z0.append(rad_i)

z0.append(param.G_t)

# Create Model
model = create_model(matrix)
print("Model Created")

# ODE Solve
z = odeint(model, z0, t)


for i in range(num_biofilms):
    print(i)
    modulo_phase_1 = np.mod(z[:,i], 2 * np.pi)
    plt.plot(t, modulo_phase_1)
    
print()

plt.show()

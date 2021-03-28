import numpy as np 
import parameters as param

def step(x):
    '''
    Step function, or Heaviside function which evaluates to 0 for any value
    less than or equal to 0 evaluates to 0. Otherwise, it is 1.

    Input: x, a number to pass into the heaviside function
    '''
    if x > 0:
        return 1
    else:
        return 0

def kuramoto(G, custom_k_0):
    '''
    Describes coupling strength for biofilms which is dependent on the amount 
    of current Glutamate in the system. K_0 represents the max coupling strength 
    k_theta represents the saturation value for oscillation coupling strength
    based on the glutamate concentration.

    Input: G, the glutamate concentration inside of a microfluidic system at
    some time, t.
    '''
    if (custom_k_0 == None):
        os = (param.K_0 * G) / (param.k_theta + G)
    else:
        os = (custom_k_0 * G) / (param.k_theta + G)
    return os 

def d_omega(G, theta):
    '''
    Expresses rate in change in the intrinsic frequency of some arbitrary biofilm 
    based on the stress level that a particular biofilm may have encountered. 
    This term increases as during the stress build up phase as G goes down. 

    Input: G, the concentration of glutamate at some time, t
    Input: theta, the phase variable of a biofilm oscillator at some time t.
    '''
    cos_out = np.cos(theta)
    denominator = 1 + (G / param.k_omega)
    out = (param.delta_w_0 / denominator) * cos_out * step(cos_out)
    return out 

def g_con(G, theta):
    '''
    Refers to glutamate consumption for a particular biofilm, representative of
    theta, the phase variable. The 1 - sin(theta) represents the increased rate
    of consumption in the case of less stress. That is, when biofilms are in 
    phase with each other, they can consume more glutamate due to the implied
    understanding that they are only in phase when there is available glutamate.

    Input: G, the glutamate concentration at a certain time, t
    Input: theta, the phase variable for a certain biofilm at a time, t
    '''
    consumption_rate = ((param.alpha * G) / (param.k_g + G))
    stress = (1 - np.sin(theta))
    return  consumption_rate * stress 

def g_add(g):
    '''
    This function is subject to change depending on the media flow through the
    microfluidic chamber.

    Input: G, the glutamate concentration at a certain time, t.
    '''
    return param.beta * g

def m_consump(r, g, custom_delta_g):
    '''
    This function defines the metabolic consumption rate of a particular biofilm
    given its size, r, and the current glutamate concentration in the microfluidic
    chamber.

    Input: r, the current size of some biofilm at time, t
    Input: g, the current glutamate concentration in microfluidic chamber
    Input: custom_delta_g, if None, use default param, but otherwise represents
    custom delta_g to use.
    '''
    if (custom_delta_g == None):
        return param.delta_g * r * g
    else:
        return custom_delta_g * r * g

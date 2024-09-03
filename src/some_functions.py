import numpy as np

def Eyring_k(G, T=300):
    k = 1.38e-23  # Boltzmann constant (J/K)
    h = 6.626e-34  # Planck's constant (J·s)
    R = 8.314  # Universal gas constant (J/(mol·K))
    
    # G is now in kJ/mol, so we directly use it with R
    return k * T / h * np.exp(-G * 1000 / (R * T))  # Convert kJ to J
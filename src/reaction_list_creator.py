import numpy as np
import pandas as pd
import joblib
from some_functions import Eyring_k

def insertion(n, k_ins, k_deins):
    return [f'CatC{n} + C2 <-> CatC{n+2}', k_ins, k_deins]

def deattaching(n, k_deattaching, k_attaching):
    return [f'CatC{n} <-> Cat_C{n}', k_deattaching, k_attaching]

def adsorption(n, k_adsorption, k_desorption):
    return [f'CatC{n} <-> Cat + C{n}', k_desorption, k_adsorption]

def deorientation(n, k_deorientation, k_orientation):
    return [f'Cat::C{n} <-> Cat + C{n}', k_deorientation, k_orientation]

def diffusion(n, k_diffusion_dict):
    return [f'C{n} <D> C{n}_out', k_diffusion_dict[f'k_diffusion_{n}'], k_diffusion_dict[f'k_diffusion_{n}']]

def adsorption_catalityc_deactiv(n, k_adsorption, k_desorption):
    return [f'CatC{n} <-> Cat* + C{n}', k_desorption, k_adsorption]

def catalityc_activ(k_activ, k_deactiv):
    return [f'Cat* <-> Cat', k_activ, k_deactiv]

def make_known_constants(barriers, T=300):
    k = {}
    e_desorption = max(barriers['c4_desorbtion'], barriers['c4_deattaching'])
    e_adsorption = max(barriers['c4_adsorption'], barriers['c4_attaching'])
    k['k_adsorption'] = Eyring_k(e_adsorption, T)
    k['k_desorption'] = Eyring_k(e_desorption, T)
    k['k_ins'] = Eyring_k(barriers['c2_insertion'], T)
    k['k_deins'] = Eyring_k(barriers['c2_deinsertion'], T)
    k['k_c2_adsorption'] = Eyring_k(barriers['c2_adsorption'], T)
    k['k_c2_desorption'] = Eyring_k(barriers['c2_desorbtion'], T)
    return k

def add_diffusion(constants, diffusion):
    d = dict()
    for n in range(2, 29, 2):
        d[f'k_diffusion_{n}'] = diffusion[n]
    constants['k_diffusion_dict'] = d
    return constants

def add_orientation(constants, orientation=1e-2, deorientation=1e2):
    constants[f'k_orientation'] = orientation
    constants[f'k_deorientation'] = deorientation
    return constants

def make_basic_reactions(n, k_adsorption, k_desorption,  
                   k_ins, k_deins, k_diffusion_dict,
                   k_c2_adsorption, k_c2_desorption):
    if n == 2:
        reactions = [adsorption(n, k_c2_adsorption, k_c2_desorption),
                     insertion(n, k_ins, k_deins), 
                    diffusion(n, k_diffusion_dict)]
    else:
        reactions = [adsorption(n, k_adsorption, k_desorption), 
                    insertion(n, k_ins, k_deins), 
                    diffusion(n, k_diffusion_dict)]
    
    
    return reactions

def make_cat_ac_reactions(n, k_adsorption, k_desorption,  
                   k_ins, k_deins, k_diffusion_dict,
                   k_c2_adsorption, k_c2_desorption):
    if n == 2:
        reactions = [adsorption_catalityc_deactiv(n, k_c2_adsorption, k_c2_desorption),
                     insertion(n, k_ins, k_deins), 
                    diffusion(n, k_diffusion_dict)]
    else:
        reactions = [adsorption_catalityc_deactiv(n, k_adsorption, k_desorption), 
                    insertion(n, k_ins, k_deins), 
                    diffusion(n, k_diffusion_dict)]
    
    return reactions


def make_reaction_list(n, k, k_ads=None, diffusion_dict=None, k_activ=None, 
                       k_deactiv=None, borders=(1e-7, 1e7)):
    
    k = k.copy()

    constants = make_known_constants(k)
    if k_ads is not None:
        constants['k_adsorption'] = k_ads
    
    for key, val in constants.items():
        constants[key] = val if val < borders[1] else borders[1]

    for key, val in constants.items():
        constants[key] = val if val > borders[0] else borders[0]

    constants = add_diffusion(constants, diffusion_dict)
    
    reaction_list = []
    if k_deactiv is None:
        for i in range(2, n+2, 2):
            reaction_list += make_basic_reactions(i, **constants)
        reaction_list = reaction_list[:-2]
    
    else:
        for i in range(2, n+2, 2):
            reaction_list += make_cat_ac_reactions(i, **constants)
        reaction_list = reaction_list[:-2]
        reaction_list += [catalityc_activ(k_activ, k_deactiv)]

    return reaction_list
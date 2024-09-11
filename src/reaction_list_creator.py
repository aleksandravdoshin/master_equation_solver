import numpy as np
import pandas as pd
import joblib
from some_functions import Eyring_k

def insertion(n, k_ins, k_deins):
    return [f'CatC{n} + C2 <-> CatC{n+2}', k_ins, k_deins]

def adsorption(n, k_adsorption, k_desorption):
    return [f'Cat + C{n} <-> CatC{n}', k_adsorption, k_desorption]

def physisorption(n, k_physisorption, k_dephysisorption, ):
    return [f'Cat + C{n} <-> Cat::C{n}', k_physisorption, k_dephysisorption]

def chemisorption(n, k_chemisorption, k_dechemisorption):
    return [f'Cat::C{n} <-> Cat_C{n}', k_chemisorption, k_dechemisorption]

def physis_chemi_sorption(n, k_phys_chem, k_chem_phys):
    return [f'Cat::C{n} <-> Cat_C{n}', k_phys_chem, k_chem_phys]

def diffusion(n, k_diffusion_dict):
    return [f'C{n} <D> C{n}_out', k_diffusion_dict[f'k_diffusion_{n}'], k_diffusion_dict[f'k_diffusion_{n}']]

def add_diffusion(n_max, constants, diffusion, ):
    d = dict()
    for n in range(2, n_max+1, 2):
        d[f'k_diffusion_{n}'] = diffusion[n]
    constants['k_diffusion_dict'] = d
    return constants

def make_constants_insertion(barriers, T=300, additional_barrier=None):
    k = {}
    c2_insertion = barriers['c2_insertion']
    c2_deinsertion = barriers['c2_deinsertion']
    if additional_barrier is not None:
        if c2_insertion > 0:
            c2_insertion += additional_barrier
        else:
            c2_insertion = additional_barrier
        if c2_deinsertion > 0:
            c2_deinsertion += additional_barrier
        else:
            c2_deinsertion = additional_barrier
    k['k_ins'] = Eyring_k(c2_insertion, T)
    k['k_deins'] = Eyring_k(c2_deinsertion, T)
    return k

def make_constants_c2_ads_des(barriers, T=300, des_barrier=None):
    k = {}
    e_c2_ads = barriers['c2_adsorption']
    e_c2_des = barriers['c2_desorbtion']
    if des_barrier is not None:
        if e_c2_ads > 0:
            e_c2_ads += des_barrier
        else:
            e_c2_ads = des_barrier
        if e_c2_des > 0:
            e_c2_des += des_barrier
        else:
            e_c2_des = des_barrier
    k['k_c2_adsorption'] = Eyring_k(e_c2_ads, T)
    k['k_c2_desorption'] = Eyring_k(e_c2_des, T)
    return k

def make_constants_c2_chem_ads_des(barriers, T=300, des_barrier=None):
    k = {}
    e_c2_ads = barriers['c2_chemisorb']
    e_c2_des = barriers['c2_dechemisorb']
    if des_barrier is not None:
        if e_c2_ads > 0:
            e_c2_ads += des_barrier
        else:
            e_c2_ads = des_barrier
        if e_c2_des > 0:
            e_c2_des += des_barrier
        else:
            e_c2_des = des_barrier
    k['k_chemisorption'] = Eyring_k(e_c2_ads, T)
    k['k_dechemisorption'] = Eyring_k(e_c2_des, T)
    return k

def make_constants_c4_ads_des(barriers, T=300, des_barrier=None):
    k = {}
    e_c4_ads = barriers['c4_adsorption']
    e_c4_des = barriers['c4_desorbtion']
    if des_barrier is not None:
        if e_c4_ads > 0:
            e_c4_ads += des_barrier
        else:
            e_c4_ads = des_barrier
        if e_c4_des > 0:
            e_c4_des += des_barrier
        else:
            e_c4_des = des_barrier
     
    k['k_adsorption'] = Eyring_k(e_c4_ads, T)
    k['k_desorption'] = Eyring_k(e_c4_des, T)
    return k

def make_constants_c4_chem_ads_des(barriers, T=300, des_barrier=None):
    k = {}
    e_c4_ads = barriers['c4_chemisorb']
    e_c4_des = barriers['c4_dechemisorb']
    if des_barrier is not None:
        if e_c4_ads > 0:
            e_c4_ads += des_barrier
        else:
            e_c4_ads = des_barrier
        if e_c4_des > 0:
            e_c4_des += des_barrier
        else:
            e_c4_des = des_barrier
    k['k_chemisorption'] = Eyring_k(e_c4_ads, T)
    k['k_dechemisorption'] = Eyring_k(e_c4_des, T)
    return k

def make_all_dft_constants(barriers, T=300, c2_ins_barrier=None, c2_barrier=None, c2_chem_barrier=None, c4_barrier=None, c4_chem_barrier=None):
    k = make_constants_insertion(barriers, T, c2_ins_barrier)
    k.update(make_constants_c2_ads_des(barriers, T, c2_barrier))
    k.update(make_constants_c2_chem_ads_des(barriers, T, c2_chem_barrier))
    k.update(make_constants_c4_ads_des(barriers, T, c4_barrier))
    k.update(make_constants_c4_chem_ads_des(barriers, T, c4_chem_barrier))
    return k

def border_constants(constants, borders):
    for key, val in constants.items():
        if key != 'k_diffusion_dict':
            constants[key] = val if val < borders[1] else borders[1]
    for key, val in constants.items():
        if key != 'k_diffusion_dict':
            constants[key] = val if val > borders[0] else borders[0]

def prepare_constants(n, k, borders=(1e-7, 1e7),diffusion_dict=None, **kwargs):
    constants = make_all_dft_constants(k, **kwargs)
    if diffusion_dict is not None:
        constants = add_diffusion(n, constants, diffusion_dict)
    border_constants(constants, borders)
    return constants

def make_reaction_list_basic(n, constants):
    reaction_list = []
    c = {'k_adsorption', 'k_desorption', 'k_ins', 
                 'k_deins', 'k_c2_adsorption', 'k_c2_desorption'}
    keys = set(constants.keys())
    keys -= c

    for i in keys:
        constants.pop(i)

    for i in range(2, n+2, 2):
        reaction_list += make_basic_reactions(i, **constants)
    reaction_list = reaction_list[:-2]

    return reaction_list

def make_basic_reactions(n, k_adsorption, k_desorption,  
                   k_ins, k_deins, k_c2_adsorption, k_c2_desorption):
    if n == 2:
        reactions = [adsorption(n, k_c2_adsorption, k_c2_desorption),
                     insertion(n, k_ins, k_deins)]
    else:
        reactions = [adsorption(n, k_adsorption, k_desorption), 
                    insertion(n, k_ins, k_deins)]
    
    return reactions
        

def make_reaction_list_with_diffusion(n, constants, diffusion_dict):
    reaction_list = make_reaction_list_basic(n, constants)
    for i in range(2, n+2, 2):
        reaction_list.append(diffusion(i, diffusion_dict))

    return reaction_list


from reaction_list_creator import *
from some_functions import *
from diffusion import Reaction_Diffusion
import re
import numpy as np

def run_model(n=28, x=(5e+08, 5.1e+01, 1.2e-05), cat=None, c2_conc=0.4, hours=1, capacity=None, diffusion_dict=None, n_steps=1e5):
    reaction_list = make_reaction_list(n, cat, k_ads=x[0], diffusion_dict=diffusion_dict,
                    k_activ=x[1], k_deactiv=x[2], borders=(1e-9, 1e9))
    r = Reaction_Diffusion(reaction_list, capacity=capacity)
    start_time = 0
    end_time = 60 * 60 * hours
    n_steps = int(n_steps)

    c2_conc = c2_conc
    result = r.solve({'Cat': 4, 'C2': c2_conc}, np.linspace(start_time, end_time, n_steps))
    return result

def calc_distribution(result):
    pat = re.compile(r'\bC(\d+)')
    columns = list(filter(lambda x: pat.match(x), result.columns.values))
    df = result[columns].iloc[-1000:]
    df = df.mean()
    df = df.reset_index(name='concentration')
    df['n'] = df['index'].str.extract(r'(\d+)')[0].astype(int)
    df = df.groupby('n').sum(numeric_only=True)
    return df


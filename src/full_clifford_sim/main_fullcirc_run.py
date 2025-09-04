import stim
import sinter
import pymatching
import numpy as np
import random
from typing import List, Callable, Tuple
import matplotlib.pyplot as plt
from dataclasses import dataclass

from full_clifford_sim.coords import *
from full_clifford_sim.full_circuit_fxns import *
from full_clifford_sim.s3_fxns import *
from full_clifford_sim.main_complied_fxns import *
from full_clifford_sim.gap_sampler import *


if __name__ == '__main__':


    gs_tasks = [
        sinter.Task(
            circuit = full_circuit(p, 
                                   dfinal=df, 
                                   ghz_size = glen,
                                   latter_rounds=l, prep=prep, ps_on_d3=1),
            json_metadata={'p': p, 'b':'Y', 'noise':'uniform', 'd2':df, 
                            'c':f'e2e-Y-{prep}-g{glen}-3ps1-stab{df}x{l}',
                            'ghz_size':3, 'latter_rounds': l},
        )
        for p in [0.001]
        for prep in ["hookinj"]
        for glen in [3]
        for l in [3]
        for df in [13]
    ]
    
    res = sinter.collect(tasks=gs_tasks, 
                    num_workers=1, 
                    decoders = ['pymatching-gap'],
                    custom_decoders=sinter_samplers(),
                     max_shots=80_000_000,  # Set a reasonable default
                     max_errors=10_000_000 ,   # Set a reasonable default
                     save_resume_filepath ="your_outfile_here.csv"
                     )
    print(res)


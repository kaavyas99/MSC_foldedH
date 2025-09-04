import stim
from full_clifford_sim.coords import *
from full_clifford_sim.full_circuit_fxns import *
from full_clifford_sim.s3_fxns import *
from full_clifford_sim.noise_model import insert_circuit_errs
import full_clifford_sim.ug_coords as sc

def full_circuit(nm: float, 
                 dfinal: int, 
                 prep: str,
                 latter_rounds: int = 3,
                 ghz_size: int = 3,
                 component_array: List = [1,1,1,1,1,1],
                 cultiv_only: bool = False,
                 ps_on_d3: int = 0,
                 neutralatom: bool = False,
                 handoff: bool= False
                 ) -> stim.Circuit:
    

    rsc = FullCircuit(dx=dfinal,
                    dy=dfinal ,
                    glen=ghz_size, 
                    basis="Y")
    
    #component_array [0:unitary prep, 1:ghz_prep+dec, 2:cbasis_check, \
    # 3:uni_grow, 4:d5_stab, 5:final_growth]

    if prep == "hookinj":
        #stab base hook inj + Y prep
        prep_circuit = rsc.cstage_circ.d3rot_hookinj()
        prep_circuit += rsc.cstage_circ.d3_rot_to_reg()
    elif prep == "optunit":
        #unitary prep of Y state
        prep_circuit  = rsc.cstage_circ.opt_rotd3_uprep()
        prep_circuit += rsc.cstage_circ.d3_modrotmeas() #only measures 3 stabs
        prep_circuit += rsc.cstage_circ.d3_rot_to_reg()
    elif prep == "unitstab":
        prep_circuit  = rsc.cstage_circ.d3_prep_circuit(to_reg=False)
        prep_circuit += rsc.cstage_circ.d3_rotmeas()
        prep_circuit += rsc.cstage_circ.d3_rot_to_reg()
    else:
        raise NotImplementedError

    rsc.qcircuit += insert_circuit_errs(prep_circuit, nm, 
        valid = component_array[0], convert_nac=neutralatom)   

    #CH check
    for _ in range(2): 

        #ghz prep
        ghz_prep = rsc.ghzcirc.prepare_ghz_state()
        rsc.qcircuit += insert_circuit_errs(ghz_prep, nm, 
            valid = component_array[1], convert_nac=neutralatom)

        #noisy check
        noisy_check = rsc.cbasis_check()
        rsc.qcircuit += insert_circuit_errs(noisy_check, nm+(nm*neutralatom), 
            valid = component_array[2], convert_nac=neutralatom)
        
        ghz_meas = rsc.ghzcirc.measure_ghz_state()
        rsc.qcircuit += insert_circuit_errs(ghz_meas, nm, 
            valid = component_array[1], convert_nac=neutralatom)
            
        
    if ps_on_d3 == 1: # we might choose to PS on Reg(3)
        stay_ps =  rsc.cstage_circ.d3reg_stabmsmt()
        rsc.qcircuit += insert_circuit_errs(stay_ps, nm,
                valid = component_array[4], convert_nac=neutralatom)  
    elif ps_on_d3 == 2: #we might choose stab msmt on Rot(3)
        print("WARNING: doing ps on Rot")
        rot_ps = rsc.cstage_circ.d3_rotmeas(onlylasttwo=True, prev=False)
        rsc.qcircuit += insert_circuit_errs(rot_ps, nm,
                valid = component_array[4], convert_nac=neutralatom)  
    
    #grow d3 to d5 using unitary encoder     
    if ps_on_d3 != 2:
        grow_d3d5 = rsc.cstage_circ.grow_3u5r()
        rsc.qcircuit += insert_circuit_errs(grow_d3d5, nm,
            valid = component_array[3], convert_nac=neutralatom)

    #do postselecion on stab msmt at d-5
    if ps_on_d3 == 0: 
        grow_ps = rsc.sc_stab_round(d_rest=5)
        grow_ps += rsc.sc_detectors(curr_only=True, 
                                    d_rest=5, 
                                    ps_round=True)
        rsc.qcircuit += insert_circuit_errs(grow_ps, nm,
            valid = component_array[4], convert_nac=neutralatom)

    if cultiv_only: #always returns at Rot(5)
        return rsc.qcircuit
        
    
    d_rest = 3 if ps_on_d3 == 2 else 5
    #reset larger code, prepare for YLi lattice surgery
    grow_final = rsc.larger_code_reset(d_rest=d_rest)

    #do first stabilizer msmt round on larger code
    grow_final += rsc.sc_stab_round()
    grow_final += rsc.sc_detectors(curr_only=True, 
                                    first_round=True, 
                                    d_rest=d_rest)

    #do latter stab msmts on larger code
    for _ in range(latter_rounds - 1):
        grow_final += rsc.sc_stab_round()
        grow_final += rsc.sc_detectors()

    #insert errors into above noisy circuits
    rsc.qcircuit += insert_circuit_errs(grow_final, nm,
        valid = component_array[5], convert_nac=neutralatom)

        
    #last perfect stab msmt round for decoding
    rsc.qcircuit += rsc.sc_stab_round()
    rsc.qcircuit += rsc.sc_detectors()

    #add logical msmt
    if not handoff:
        rsc.qcircuit += rsc.logYMeas()
        print(len(rsc.qcircuit.shortest_graphlike_error()))
    else:
        #add logical msmt (stim 1.15)
        z_targs = [2*i for i in range(0,rsc.dx)]
        x_targs = [2*rsc.dx*i for i in range(0,rsc.dy)]
        rsc.qcircuit.append("OBSERVABLE_INCLUDE", arg=0, 
                            targets=[stim.target_z(i) for i in z_targs])
        rsc.qcircuit.append("OBSERVABLE_INCLUDE", arg=1, 
                            targets=[stim.target_x(i) for i in x_targs])

    return rsc.qcircuit
    
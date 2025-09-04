from dataclasses import dataclass
import stim
from full_clifford_sim.coords import *

@dataclass
class GHZstate:
    dx: int
    dy: int
    glen: int


    def layout_ghz_state(self) -> stim.Circuit:
        "Lay out a GHZ state"

        ghz_circ = stim.Circuit()
        self.sc_offset =  2*(self.dx+1)*(self.dy+1)

        ghz_qubs = np.array(range(self.glen))  + self.sc_offset

        #lay out ghz state
        for qidx, q in enumerate(ghz_qubs):
            q_x = qidx
            q_y = self.dy + 2 
            ghz_circ.append("QUBIT_COORDS",  
                q, (q_x, q_y,0))
        #add flag qubit coordinates
        if self.glen > 3:
            ghz_circ.append("QUBIT_COORDS",  
                    q+1, (q_x+1, q_y,0))
            

        self.ghz_qubs = ghz_qubs
        return ghz_circ
    
    def prepare_ghz_state(self) -> stim.Circuit:
        "Entangle a GHZ state"

        ghz_circ = stim.Circuit()
        mid_q = self.glen//2  + self.sc_offset

        if self.glen == 5: 
            gate_array = [[mid_q, mid_q + 1], 
                          [mid_q, mid_q - 1, mid_q + 1, mid_q+ 2],
                          [mid_q-1, mid_q - 2] ]
                           
        elif self.glen == 3:
            gate_array = [[mid_q, mid_q + 1], 
                          [mid_q, mid_q - 1 ]]
                          
        else:
            raise NotImplementedError 
            
        self.prep_gate_array = gate_array
            
        ghz_circ.append("R", self.ghz_qubs)
        ghz_circ.append("TICK")
        ghz_circ.append("H", mid_q)
        ghz_circ.append("TICK")
    
        #do CNOT gates for ghz state
        for qidx in range(self.glen//2 + 1):
            filtered_gates =  gate_array[qidx]
            ghz_circ.append("CNOT", filtered_gates)
            ghz_circ.append("I", list(set(self.ghz_qubs) - set(filtered_gates)))
            ghz_circ.append("TICK")
                
    
        if self.glen == 5:
            # add flag
            last_qub = self.ghz_qubs[-1]
            first_qub = self.ghz_qubs[0]
            ghz_circ.append("CNOT", [last_qub, last_qub+1])
            ghz_circ.append("I", self.ghz_qubs[:-1])
            ghz_circ.append("TICK")
    
            ghz_circ.append("CNOT", [first_qub, last_qub+1])
            ghz_circ.append("I", self.ghz_qubs[1:])
            ghz_circ.append("M", last_qub+1)
            ghz_circ.append("TICK")
    
            ghz_circ.append("DETECTOR", stim.target_rec(-1))

        return ghz_circ
    

    def measure_ghz_state(self) -> stim.Circuit:
        "Measure and detect a GHZ state"

        ghz_circ = stim.Circuit()
        mid_q = self.glen//2  + self.sc_offset

        #do CNOT gates for ghz state
        for qidx in range(self.glen//2,-1,-1):
            filtered_gates =  self.prep_gate_array[qidx]
            ghz_circ.append("CNOT", filtered_gates)
            ghz_circ.append("I", list(set(self.ghz_qubs) - set(filtered_gates)))
            ghz_circ.append("TICK")

        ghz_circ.append("H", mid_q)
        ghz_circ.append("TICK")

        ghz_circ.append("M", self.ghz_qubs)
        for idx in range(self.glen):
            ghz_circ.append("DETECTOR", stim.target_rec(-idx-1))

        ghz_circ.append("SHIFT_COORDS", arg=(0,0,1))
        ghz_circ.append("TICK")
            
        return ghz_circ
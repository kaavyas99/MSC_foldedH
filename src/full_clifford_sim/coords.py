from dataclasses import dataclass
import stim
import numpy as np
from typing import Tuple, List

@dataclass
class RotSurfCodeCoords:
    "Qubit positions and gate pairs for one large surface code"
    dx: int
    dy: int

    def is_valid_ancilla(self,a, d_rest : int = None) -> bool:
        "Origin is top left, with top left two body stabilizer on top"
        if a < 0: return False
        
        a = a //2

        if d_rest is None:
            dx, dy = self.dx, self.dy
        else:
            dx, dy = d_rest, d_rest

        xind = a %  (self.dx + 1) 
        yind = a // (self.dx + 1)

        if yind == 0: # if on top
            if xind % 2 == 1: return False
        if xind == 0: # if on left
            if yind % 2 == 0: return False
        if yind == dy: # if on  bottom
            if xind % 2 == 0: return False
        if xind == dx: # if on right
            if yind % 2 == 1: return False
        if xind > dx: return False #outside scope of code
        if yind > dy: return False

        return True

    
    def is_Xstab(self, a) -> bool:
        "returns True for X stabilizers"
        a = (a % self.no1SCQubits)//2
        posParity = a % (self.dx +1 ) + a // (self.dx+1)
        if posParity % 2 == 1 : return False
        return True
    

    def q_parity(self,q) -> bool:
        "Parity helps determine gate ordering"
        q = (q % self.no1SCQubits)/2
        return bool(((q % self.dx) + (q // self.dx))%2)


    def nw_ancilla(self,q) -> int:
        "Get ancilla SW to data qubit"
        qsc = q // self.no1SCQubits
        q = (q % self.no1SCQubits)/2

        scA = (self.dx +1)  * ( q // self.dx )
        scA += q % self.dx 

        return int(qsc * self.no1SCQubits + 2*scA +1)


    def __post_init__(self):
        "Initialize qubit and gate arrays"
        dx = self.dx
        dy = self.dy

        self.no1SCQubits =  2*(dx+1)*(dy+1)

        self.dataQubits = 2*np.arange(dx*dy)
        
        roughAncillas = 2*np.arange((dx+1)*(dy+1))+ 1 
        
        self.ancillaQubits = np.array([a for a in roughAncillas 
                                       if self.is_valid_ancilla(a)])
        
        self.XancillaQubits = np.array([a for a in self.ancillaQubits
                                       if self.is_Xstab(a)])
        
        self.ZancillaQubits = np.array([a for a in self.ancillaQubits
                                       if not self.is_Xstab(a)])

        self.allQubits = np.array(list(self.ancillaQubits) 
                                  + list(self.dataQubits))

        #use a CNOT-only gateset
        stepGates = [[],[],[],[]]

        for q in self.dataQubits:
            nwA = self.nw_ancilla(q)
            swA = nwA + 2*(self.dx+1)
            neA = nwA + 2
            seA = swA + 2

            parityOrder = [[swA,seA,nwA,neA],[swA,nwA,seA,neA]]
            qOrder = self.q_parity(q)

            for step in range(4):
                stepA = parityOrder[qOrder][step]
                if stepA in self.ancillaQubits:
                    if self.is_Xstab(stepA):
                        stepGates[step] += [stepA,q]
                    elif not self.is_Xstab(stepA):
                        stepGates[step] += [q, stepA]


        self.stepGates = np.array(stepGates)

        idleQubits = [[],[],[],[]]
        for step in range(4):
            idleQubits[step] = list(set(self.allQubits) 
                                    - (set(stepGates[step])))
        self.idleQubits = idleQubits

    
    def layout_coords(self) -> stim.Circuit:
        "Initialize a stim circuit with the qubits in SC config"

        scCirc = stim.Circuit()

        for dq in self.dataQubits:
            qc = dq/2
            qc_x = qc % self.dx 
            qc_y = qc //self.dx 
            scCirc.append("QUBIT_COORDS",  
                dq, (qc_x, qc_y,0))

        for aq in self.ancillaQubits:
            scCirc.append(
               "QUBIT_COORDS", 
                aq, 
                self.ancilla_coords(aq))

        return scCirc

    def ancilla_coords(self, aq: int) -> Tuple:
        "Returns a tuple of ancilla coordinates"
        qc = aq //2
        qc_x = qc % (self.dx+1) - 0.5
        qc_y = qc // (self.dx+1) - 0.5
        return (qc_x, qc_y,0)
    
    
    
    def filter_step_gates(self, d_rest: int) -> List:
        "take stepGates and return stepGates for a smaller d_rest x drest SC"

        filtered_step_gates =[]

        for this_step_gates in self.stepGates:
            edited_sg = []
            for i in range(len(this_step_gates)// 2):
                q1 = this_step_gates[2*i]
                q2 = this_step_gates[2*i +1]
                dq = q1 if q1 % 2 == 0 else q2
                aq = q1 if q1 % 2 != 0 else q2

                dq_x = (dq / 2) % self.dx 
                dq_y = (dq / 2) //self.dx 

                if (dq_x < d_rest and dq_y < d_rest):
                    if self.is_valid_ancilla(aq, d_rest=d_rest):
                        edited_sg.extend([q1, q2])

            filtered_step_gates.append(edited_sg)

        return filtered_step_gates    



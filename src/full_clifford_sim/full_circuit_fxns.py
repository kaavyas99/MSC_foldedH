from dataclasses import dataclass
import stim
from full_clifford_sim.coords import *
from full_clifford_sim.ghz_fxns import *
from full_clifford_sim.s3_fxns import *

@dataclass
class FullCircuit:
    dx: int
    dy: int
    glen: int
    basis: str
    smallsc: bool

    def sc_stab_round(self, d_rest: int = None) -> stim.Circuit:
        "Returns gates and ancilla msmts for any stab rd"

        if d_rest is None:
            stepGates = self.cbase.stepGates
        else:
            stepGates = self.cbase.filter_step_gates(d_rest)

        fround = stim.Circuit()
        fround.append("R", self.cbase.ancillaQubits)
        fround.append("TICK")

        fround.append("H", self.cbase.XancillaQubits)
        fround.append("TICK")

        for step in range(4):
            fround.append("CNOT", stepGates[step])
            fround.append("TICK")

        fround.append("H", self.cbase.XancillaQubits)
        fround.append("TICK")

        #measure ancillas, idling noise on datas
        fround.append("M", self.cbase.ancillaQubits)
        fround.append("TICK")
        
        return fround
    

    
    def sc_detectors(self, 
                     d_rest: int = None,
                     curr_only: bool = False,
                     first_round: bool = False,
                     ps_round: bool = False) -> stim.Circuit:
        "Defines a regular detector set for the input SC"

        detCircuit = stim.Circuit()
        if d_rest is None: d_rest = self.dx

        no1scA = self.no1SCAnc

        for k in range(no1scA):
            aq = self.cbase.ancillaQubits[k]
            det_loc_arg =  (*self.cbase.ancilla_coords(aq), ps_round)
            if not first_round and self.cbase.is_valid_ancilla(aq, d_rest=d_rest):
                if not curr_only:
                    detCircuit.append("DETECTOR", 
                                  [stim.target_rec(-no1scA+k),
                                   stim.target_rec(-len(self.cbase.ancillaQubits)-(no1scA)+k)],
                                   arg = det_loc_arg) 
                else: #current only
                    if not first_round:
                        detCircuit.append("DETECTOR", 
                                    [stim.target_rec(-no1scA+k)],
                                    arg = det_loc_arg) 
                        
            if first_round and curr_only: #custom detectors for transition round
                qc = aq //2
                qc_x = qc % (self.dx+1) - 0.5
                qc_y = qc // (self.dx+1) - 0.5
                if (qc_y > d_rest and qc_y > qc_x and self.cbase.is_Xstab(aq)):
                        detCircuit.append("DETECTOR", 
                            [stim.target_rec(-no1scA+k)],
                            arg = det_loc_arg) 
                elif (qc_x > d_rest and qc_x > qc_y  and not self.cbase.is_Xstab(aq)):
                        detCircuit.append("DETECTOR", 
                            [stim.target_rec(-no1scA+k)],
                            arg = det_loc_arg) 
                elif self.cbase.is_valid_ancilla(aq, d_rest=d_rest):
                    detCircuit.append("DETECTOR", 
                            [stim.target_rec(-no1scA+k)],
                            arg = det_loc_arg) 

        detCircuit.append("SHIFT_COORDS", arg=(0,0,1))

        return detCircuit

    

    def larger_code_reset(self, d_rest: int = 5) -> stim.Circuit:
        "initialize larger code in Ying Li- like state"
        r"""
        ------ --
        | psi |0|
        _______ |
        | +    \|
        _________
        
        """
        yli_circ = stim.Circuit()

        outside_qubs = []
        plus_qubs = []

        #first reset all data qubits outside OG surface code
        for q in self.cbase.dataQubits:
            qc = q/2
            qc_x = qc % self.dx 
            qc_y = qc //self.dx 

            if (qc_x >= d_rest or qc_y >= d_rest):
                outside_qubs.append(q)
                #now put lower qubits in |+> state
                if (qc_y >= d_rest and qc_y > qc_x):
                    plus_qubs.append(q)

        yli_circ.append("R", outside_qubs)
        yli_circ.append("TICK")
        yli_circ.append("H", plus_qubs)
        yli_circ.append("TICK")

        return yli_circ

    
    def logYMeas(self) -> stim.Circuit:
        "returns logical Y msmt and obs"

        measCircuit = stim.Circuit()
        yMeas = []
        
        measCircuit.append("MY", 0)
        measCircuit.append("MZ", [2*i for i in range(1,self.dx)])
        measCircuit.append("MX", [2*self.dx*i for i in range(1,self.dy)])
        measCircuit.append("TICK")

        measlen = self.dx + self.dy - 1

        yMeas = [ stim.target_rec(-measlen + k)
                 for k in range(measlen)] 

        measCircuit.append("OBSERVABLE_INCLUDE", yMeas, 0) 
        return measCircuit
    
    
    def cbasis_check(self) -> stim.Circuit:
        "does a transvesral check using the GHZ state"

        if self.basis not in ["Y"]:
            raise ValueError(f"The basis {self.basis} is not supported")
        
        check_circ = stim.Circuit()

        gqs = self.ghzcirc.ghz_qubs
        mqs = self.cstage_circ.mgrid_qubs
        aqs = self.cstage_circ.agrid_qubs
        
        if self.glen==5:
            cy_step1 = [gqs[0], mqs[0], gqs[2], mqs[4], gqs[4], mqs[8]]
            check_circ.append("CY", cy_step1)
            check_circ.append("TICK")

            cy_step2 = [gqs[1], mqs[3], gqs[3], mqs[6], gqs[4], mqs[7]]
            check_circ.append("CY", cy_step2)
            check_circ.append("TICK")

            cy_step3 = [gqs[1], mqs[1], gqs[3], mqs[2], gqs[4], mqs[5]]
            check_circ.append("CY", cy_step3)
            check_circ.append("TICK")
            
        elif self.glen==3:
            cy_step1 = [gqs[0], mqs[0], gqs[1], mqs[5], gqs[2], mqs[8]]
            check_circ.append("CY", cy_step1)
            check_circ.append("TICK")

            cy_step2 = [gqs[0], mqs[3], gqs[1], mqs[6], gqs[2], mqs[7]]
            check_circ.append("CY", cy_step2)
            check_circ.append("TICK")

            cy_step3 = [gqs[0], mqs[1], gqs[1], mqs[2], gqs[2], mqs[4]]
            check_circ.append("CY", cy_step3)
            check_circ.append("TICK")

        return check_circ 
    


    def __post_init__(self):

        self.no1SCAnc = self.dx * self.dy - 1

        self.cbase = RotSurfCodeCoords(self.dx,self.dy)
        self.qcircuit = self.cbase.layout_coords()

        self.ghzcirc = GHZstate(self.dx, self.dy, self.glen)
        self.qcircuit += self.ghzcirc.layout_ghz_state()

        self.cstage_circ = CultStage(self.dx,self.dy, "Y", self.smallsc)
        
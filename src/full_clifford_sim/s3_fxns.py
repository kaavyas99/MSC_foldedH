from dataclasses import dataclass
import stim
from full_clifford_sim.coords import *
import json
from itertools import chain

def flatten(xss):
    return [x for xs in xss for x in xs]

#functions that use sv coords
    # Converts the JSON indices to stim indices.
    # keys: statevector indices. entries: stim coordinates in (x,y)

json_to_coord = {'0':(0,0),
                        '1':(2,0),
                        '2':(4,0),
                        '3':(1,1),
                        '4':(3,1),
                        '5':(0,2),
                        '6':(2,2),
                        '7':(4,2),
                        '8':(1,3),
                        '9':(3,3),
                        '10':(0,4),
                        '11':(2,4),
                        '12':(4,4)}

# convert to stim coordinates
def coord_to_stim(coord:tuple[int,int],df:int)->int:
    x,y = coord
    return (y * df + x)*2



@dataclass
class CultStage:
    dx: int
    dy: int
    basis: str
    smallsc: bool

    #support ps on rot(3)
    def qub_topleft(self,q:int) -> int:
            if not self.smallsc:
                return q - 2*self.dx - 2
            else:
                assert q % 2 == 1, "q not ancilla"
                return ((q // 2 // (self.dx+1))-1) * 2* self.dx + 2*((q//2) % (self.dx+1)) -2
        
    def qub_bottomleft(self,q:int) -> int:
        if not self.smallsc:
            return q + 2*self.dx - 2
        else:
            assert q % 2 == 1, "q not ancilla"
            return (q // 2 // (self.dx+1)) * 2* self.dx + 2*((q//2) % (self.dx+1)) -2
    
    def qub_topright(self,q:int) -> int:
        if not self.smallsc:
            return q - 2*self.dx + 2
        else:
            assert q % 2 == 1, "q not ancilla"
            return ((q // 2 // (self.dx+1))-1) * 2* self.dx + 2*((q//2) % (self.dx+1))

    def qub_bottomright(self,q:int) -> int:
        if not self.smallsc:
            return q + 2*self.dx + 2
        else:
            assert q % 2 == 1, "q not ancilla"
            return (q // 2 // (self.dx+1)) * 2* self.dx + 2*((q//2) % (self.dx+1))
        
    
    def qub_above(self, q:int) -> int:
            return q - 2*self.dx
        
    def qub_toright(self, q:int) -> int:
        return q + 2
    
    def qub_below(self, q:int) -> int:
            return q + 2*self.dx
        
    def qub_toleft(self, q:int) -> int:
        return q - 2


    def opt_rotd3_uprep(self) -> stim.Circuit:
        "use mqt package to get a gate optimal uprep circ"
        uprep_circ = stim.Circuit()

        mainqubs = [ row * 4 * self.dx + 4 * col 
                    for row in range(3) for col in range(3)]
        
        altqubs = [4 * self.dx * (row) + 2*self.dx + 2 + 4 * col 
                   for row in range(2) for col in range(2)]

        uprep_circ.append("R", mainqubs + altqubs)
        uprep_circ.append("TICK")


        uprep_circ.append("H_"+self.basis+"Z", mainqubs[0])
        uprep_circ.append("H", [mainqubs[1], mainqubs[2], mainqubs[6], mainqubs[8]])
        uprep_circ.append("H", [altqubs[0], altqubs[3]])
        uprep_circ.append("TICK")

        #first round gates 
        uprep_circ.append("CNOT", [mainqubs[2], mainqubs[0], mainqubs[8], mainqubs[5]])
        uprep_circ.append("TICK")

        #2nd round gates 
        uprep_circ.append("CNOT", [mainqubs[0], mainqubs[7], mainqubs[2], mainqubs[4]])
        uprep_circ.append("TICK")

        #3rd round gates 
        uprep_circ.append("CNOT", [mainqubs[0], mainqubs[3], mainqubs[5], mainqubs[4], mainqubs[6], mainqubs[7]])
        uprep_circ.append("TICK")

        #4th round gates 
        uprep_circ.append("CNOT", [mainqubs[1], mainqubs[2], mainqubs[4], mainqubs[6]])
        uprep_circ.append("TICK")

        return uprep_circ

    
    
    def d3_prep_circuit(self, to_reg: bool= True) -> stim.Circuit:
        "unitary encoder for unrotated or rotated surface code"

        # diagram notation: qub_index-partition
        r"""

        Index with_anc == False:
        00-? ───── 02-? ───── 04-?
        │ \      /    \     / │ 
        │   02-?       04-?   │    // 2(d+1) + 1
        │  /     \    /     \ │
        00-? ───── 02-? ───── 04-? // 4d
        │ \      /   \      / │ 
        │   02-?       04-?   │    // 4(d+1) + 1
        │ /      \   /      \ │
        00-? ───── 02-? ───── 04-? // 8d
        """
        uprep_circ = stim.Circuit()

        mainqubs = list(self.mgrid_qubs)
        altqubs = list(self.agrid_qubs)
        
        # we need to use position indexing 
        # in the above arrays to get their SC coords

        uprep_circ.append("R", mainqubs + altqubs)
        uprep_circ.append("TICK")


        uprep_circ.append("H_"+self.basis+"Z", mainqubs[4])
        uprep_circ.append("H", [mainqubs[2], mainqubs[3], mainqubs[6], mainqubs[8]])
        uprep_circ.append("H", [altqubs[0], altqubs[3]])
        uprep_circ.append("TICK")

        #first round gates (and flag prep if valid)
        uprep_circ.append("CNOT", [mainqubs[2], mainqubs[1], mainqubs[3], mainqubs[0], mainqubs[4], mainqubs[7], mainqubs[8], mainqubs[5]])
        uprep_circ.append("TICK")


        #second round gates
        uprep_circ.append("CNOT", [mainqubs[3], mainqubs[4], mainqubs[7], mainqubs[1]])
        uprep_circ.append("TICK")

        #third round gates
        uprep_circ.append("CNOT", [mainqubs[3], mainqubs[1], mainqubs[6], mainqubs[7], mainqubs[8], mainqubs[4]])
        uprep_circ.append("TICK")

        #fourth round rot preop and fifth round gates - ur to rot
        uprep_circ.append("CNOT", [mainqubs[8], mainqubs[7]])
        if to_reg:
            uprep_circ += self.d3_rot_to_reg()
        else:
            uprep_circ.append("TICK")

        return uprep_circ
    
    
    def d3_rot_to_reg(self) -> stim.Circuit:
        "rotated to regular d3 SC"

        gcirc = stim.Circuit()
        altqubs = self.agrid_qubs
        mainqubs = self.mgrid_qubs

        gcirc.append("R", self.agrid_qubs)
        gcirc.append("TICK")

        gcirc.append("H", self.agrid_qubs[[0,3]])
        gcirc.append("TICK")

        gcirc.append("CNOT", [altqubs[0], mainqubs[0], mainqubs[1], altqubs[1], mainqubs[3], altqubs[2], altqubs[3], mainqubs[4]])
        gcirc.append("TICK")

        gcirc.append("CNOT", [altqubs[0], mainqubs[1], altqubs[3], mainqubs[5], mainqubs[4], altqubs[1], mainqubs[6], altqubs[2]])
        gcirc.append("TICK")

        return gcirc
    

    def d3rot_hookinj(self, skip_gauge_fix: bool = False) -> stim.Circuit:
        "hook injection followed by ugrwoth to reg(3)"

        if skip_gauge_fix:
            print('WARNING: skip_gauge_fix==True')

        hi_circ = stim.Circuit()
        bdry_ancillas = self.bdry_ancillas

        hi_circ.append("R", self.mgrid_qubs)
        hi_circ.append("R", self.agrid_qubs)
        hi_circ.append("R", bdry_ancillas)
        hi_circ.append("TICK")

        #init with + in diag
        hi_circ.append("H", self.mgrid_qubs[[1,2,4,5,6,7,8]])
        hi_circ.append("H", bdry_ancillas[[0,3]])
        hi_circ.append("H", self.agrid_qubs[[0,3]])
        hi_circ.append("TICK")

        #get gate sequences
        
        step1gatesCX = flatten([[q,self.qub_topleft(q)] for q in self.agrid_qubs[[0,3]]])
        step1gatesCX += flatten([[self.qub_topleft(q),q] for q in self.agrid_qubs[[1,2]]])
        if not self.smallsc:
            step1gatesCX += [bdry_ancillas[2] - 2*self.dx, 
                             bdry_ancillas[2], bdry_ancillas[3], bdry_ancillas[3]-2]
        else:
            step1gatesCX += [self.qub_topleft(bdry_ancillas[2]),
                              bdry_ancillas[2], bdry_ancillas[3], self.qub_topleft(bdry_ancillas[3])]
        hi_circ.append("CX", step1gatesCX)
        hi_circ.append("TICK")

        step2gatesCX = flatten([[q,self.qub_topright(q)] for q in self.agrid_qubs[[0,3]]])
        step2gatesCX += flatten([[self.qub_topright(q),q] for q in self.agrid_qubs[[1, 2]]])  ## add anc2 here
        if not self.smallsc:
            step2gatesCX += [bdry_ancillas[3], bdry_ancillas[3] + 2, 
                             bdry_ancillas[2] + 2*self.dx, bdry_ancillas[2]]
        else:
            step2gatesCX += [bdry_ancillas[3], self.qub_topright(bdry_ancillas[3]),
                              self.qub_bottomleft(bdry_ancillas[2]), bdry_ancillas[2]]
        hi_circ.append("CX", step2gatesCX)
        hi_circ.append("TICK")
        hi_circ.append("S", 2*self.dx+6 + int(self.smallsc))
        hi_circ.append("TICK")

        step3gatesCX = flatten([[q,self.qub_bottomleft(q)] for q in self.agrid_qubs[[0,3]]])
        step3gatesCX += flatten([[self.qub_bottomleft(q),q] for q in self.agrid_qubs[[1,2]]])  ## change direction
        if not self.smallsc:
            step3gatesCX += [bdry_ancillas[0], bdry_ancillas[0] - 2, 
                             bdry_ancillas[1] - 2*self.dx, bdry_ancillas[1]]
        else:
            step3gatesCX += [bdry_ancillas[0], self.qub_bottomleft(bdry_ancillas[0]), 
                             self.qub_topright(bdry_ancillas[1]), bdry_ancillas[1]]
        hi_circ.append("CX", step3gatesCX)
        hi_circ.append("TICK")


        step4gatesCX = flatten([[q,self.qub_bottomright(q)] for q in self.agrid_qubs[[0,3]]])
        step4gatesCX += flatten([[self.qub_bottomright(q),q] for q in self.agrid_qubs[[1,2]]])
        if not self.smallsc:
            step4gatesCX += [bdry_ancillas[0], bdry_ancillas[0]+ 2, 
                             bdry_ancillas[1] + 2*self.dx, bdry_ancillas[1]]
        else:
            step4gatesCX += [bdry_ancillas[0], self.qub_bottomright(bdry_ancillas[0]), 
                             self.qub_bottomright(bdry_ancillas[1]), bdry_ancillas[1]]
        hi_circ.append("CX", step4gatesCX)
        hi_circ.append("TICK")

        #measure out ancillas
        #hi_circ.append("I", self.mgrid_qubs)
        hi_circ.append("H", bdry_ancillas[[0,3]])
        hi_circ.append("H", self.agrid_qubs[[0,3]])
        hi_circ.append("TICK")

        hi_circ.append("M", self.agrid_qubs)
        hi_circ.append("M", bdry_ancillas)
        for meas in [-1,-3,-4,-5]: 
            hi_circ.append("DETECTOR", stim.target_rec(meas))
        hi_circ.append("TICK")

        hi_circ += self.d3_rotmeas(prev=True)

        

        #fast feedback and reset
        if not skip_gauge_fix:
            hi_circ += stim.Circuit(f"CX rec[-2] {(4*self.dx+4) * (2 - int(self.smallsc))}")
            hi_circ += stim.Circuit(f"CX rec[-6] {(4*self.dx) * (2 - int(self.smallsc))}")
            hi_circ += stim.Circuit(f"CX rec[-7] {(2*self.dx+4) * (2 - int(self.smallsc))} rec[-7] {(4*self.dx+4) * (2 - int(self.smallsc))}") ## This is the right gauge op, PK250620
            hi_circ += stim.Circuit(f"CZ rec[-8] 0")
            hi_circ.append("TICK")

        return hi_circ
    
    

    def d3_rotmeas(self, prev=False, onlylasttwo:bool=False) -> stim.Circuit:
        "do stab meas at rot 3 and then grow to reg 3"

        hi_circ = stim.Circuit()
        bdry_ancillas = self.bdry_ancillas
        
        step1gatesCX = flatten([[q,self.qub_topleft(q)] for q in self.agrid_qubs[[0,3]]])
        step1gatesCX += flatten([[self.qub_topleft(q),q] for q in self.agrid_qubs[[1,2]]])
        if not self.smallsc:
            step1gatesCX += [bdry_ancillas[2] - 2*self.dx, bdry_ancillas[2], bdry_ancillas[3], bdry_ancillas[3]-2]
        else:
            step1gatesCX += [self.qub_topleft(bdry_ancillas[2]), bdry_ancillas[2], 
                             bdry_ancillas[3], self.qub_topleft(bdry_ancillas[3])]


        step4gatesCX = flatten([[q,self.qub_bottomright(q)] for q in self.agrid_qubs[[0,3]]])
        step4gatesCX += flatten([[self.qub_bottomright(q),q] for q in self.agrid_qubs[[1,2]]])
        if not self.smallsc:
            step4gatesCX += [bdry_ancillas[0], bdry_ancillas[0]+ 2, bdry_ancillas[1] + 2*self.dx, bdry_ancillas[1]]
        else:
            step4gatesCX += [bdry_ancillas[0], self.qub_bottomright(bdry_ancillas[0]), 
                             self.qub_bottomright(bdry_ancillas[1]), bdry_ancillas[1]]

        ##hi_circ.append("I", self.mgrid_qubs) can be done simul
        hi_circ.append("R", self.agrid_qubs)
        hi_circ.append("R", bdry_ancillas)
        hi_circ.append("TICK")


        hi_circ.append("H", bdry_ancillas[[0,3]])
        hi_circ.append("H", self.agrid_qubs[[0,3]])
        hi_circ.append("TICK")

        hi_circ.append("CX", step1gatesCX)
        hi_circ.append("TICK")

        step2gatesCX = flatten([[q,self.qub_topright(q)] for q in self.agrid_qubs[[0,3]]])
        step2gatesCX += flatten([[self.qub_bottomleft(q),q] for q in self.agrid_qubs[[1,2]]])
        if not self.smallsc:
            step2gatesCX += [bdry_ancillas[3], bdry_ancillas[3] + 2, bdry_ancillas[2] + 2*self.dx, bdry_ancillas[2]]
        else:
            step2gatesCX += [bdry_ancillas[3], self.qub_topright(bdry_ancillas[3]), 
                             self.qub_bottomleft(bdry_ancillas[2]), bdry_ancillas[2]]
        

        hi_circ.append("CX", step2gatesCX)
        hi_circ.append("TICK")

        if onlylasttwo: 
            hi_circ = stim.Circuit()
            hi_circ.append("R", bdry_ancillas)
            hi_circ.append("TICK")
            hi_circ.append("H", [bdry_ancillas[0], bdry_ancillas[3]])
            hi_circ.append("TICK")

        step3gatesCX = flatten([[q,self.qub_bottomleft(q)] for q in self.agrid_qubs[[0,3]]])
        step3gatesCX += flatten([[self.qub_topright(q),q] for q in self.agrid_qubs[[1,2]]])
        if not self.smallsc:
            step3gatesCX += [bdry_ancillas[0], bdry_ancillas[0] - 2, bdry_ancillas[1] - 2*self.dx, bdry_ancillas[1]]
        else:
            step3gatesCX += [bdry_ancillas[0], self.qub_bottomleft(bdry_ancillas[0]), 
                             self.qub_topright(bdry_ancillas[1]), bdry_ancillas[1]]
        

        hi_circ.append("CX", step3gatesCX)
        hi_circ.append("TICK")

        hi_circ.append("CX", step4gatesCX)
        hi_circ.append("TICK")

        hi_circ.append("H", bdry_ancillas[[0,3]])
        hi_circ.append("H", self.agrid_qubs[[0,3]])
        hi_circ.append("TICK")

        hi_circ.append("M", self.agrid_qubs)
        hi_circ.append("M", bdry_ancillas)
        measqubs = list(self.agrid_qubs) + list(bdry_ancillas)
        for meas in range(1,9): 
            locarg = (measqubs[meas-1] // (2*self.dx), (measqubs[meas-1] % (2*self.dx)) //2,  0,1)
            if prev:
                hi_circ.append("DETECTOR", [stim.target_rec(-meas), stim.target_rec(-meas-8)],
                               arg = locarg)
            else:
               hi_circ.append("DETECTOR", [stim.target_rec(-meas)], arg = locarg) 
        hi_circ.append("SHIFT_COORDS", arg=(0,0,1))
        hi_circ.append("TICK")

        return hi_circ

    def d3_modrotmeas(self) -> stim.Circuit:
        "do SOME stab meas at rot 3, dictated by uprep circ . turns out only need 3 for optimal uprep"
        assert not self.smallsc, "Optimized unitary encoding not supported for ps at Rot(3)"

        hi_circ = stim.Circuit()

        used_x_anc = self.agrid_qubs[[0]]
        used_z_anc = self.agrid_qubs[[1,2]]
        used_anc = self.agrid_qubs[:3]
        
        step1gatesCX = flatten([[q,self.qub_topleft(q)] for q in used_x_anc])
        step1gatesCX += flatten([[self.qub_topleft(q),q] for q in used_z_anc])

        step4gatesCX = flatten([[q,self.qub_bottomright(q)] for q in used_x_anc])
        step4gatesCX += flatten([[self.qub_bottomright(q),q] for q in used_z_anc])

        hi_circ.append("R", used_anc)
        hi_circ.append("TICK")

        hi_circ.append("H", used_x_anc)
        hi_circ.append("TICK")

        hi_circ.append("CX", step1gatesCX)
        hi_circ.append("TICK")

        #stab meas round 2
        step2gatesCX = flatten([[q,self.qub_topright(q)] for q in used_x_anc])
        step2gatesCX += flatten([[self.qub_bottomleft(q),q] for q in used_z_anc])
        hi_circ.append("CX", step2gatesCX)
        hi_circ.append("TICK")

        step3gatesCX = flatten([[q,self.qub_bottomleft(q)] for q in used_x_anc])
        step3gatesCX += flatten([[self.qub_topright(q),q] for q in used_z_anc])
        hi_circ.append("CX", step3gatesCX)
        hi_circ.append("TICK")

        hi_circ.append("CX", step4gatesCX)
        hi_circ.append("TICK")

        hi_circ.append("H", used_x_anc)
        hi_circ.append("TICK")

        hi_circ.append("M", used_anc)
        measqubs = list(used_anc)
        for meas in range(1,4): 
            locarg = (measqubs[meas-1] // (2*self.dx), (measqubs[meas-1] % (2*self.dx)) //2,  0,1)
            hi_circ.append("DETECTOR", [stim.target_rec(-meas)], arg = locarg) 
        hi_circ.append("SHIFT_COORDS", arg=(0,0,1))
        hi_circ.append("TICK")

        return hi_circ
    

    def d3reg_stabmsmt(self) -> stim.Circuit:
        "do stab msmt on d3 unrotated -> rn configured to be noiseless"
        
        smcirc = self.grow_3u5r() #unitary encoder has first two steps

        #next two steps
        smcirc.append("CNOT", flatten([[self.qub_toright(q), q] for q in self.mgrid_qubs[[0,3,6,1,4,7]] ]))
        smcirc.append("CNOT", flatten([[q, self.qub_toleft(q)] for q in self.agrid_qubs ]))
        smcirc.append("TICK")

        smcirc.append("CNOT", flatten([[q, self.qub_above(q)] for q in self.mgrid_qubs[3:] ]))
        smcirc.append("CNOT", flatten([[self.qub_above(q), q] for q in self.agrid_qubs ]))
        smcirc.append("TICK")

        #measure out ancillas
        for row in range(3):
            firstanc = 4*self.dx*row + 2
            smcirc.append("H", [firstanc, firstanc+4])
        smcirc.append("TICK")
        
        for row in range(3):
            firstanc = 4*self.dx*row + 2
            smcirc.append("M", [firstanc, firstanc+4])
            if row!= 2:
                smcirc.append("M", [firstanc - 2 + 2*self.dx, firstanc+2 + 2*self.dx, firstanc+6+ 2*self.dx])
        smcirc.append("TICK")

        for k in range(12):
            smcirc.append("DETECTOR", stim.target_rec(-1-k))

        return smcirc
    
    def reg_3_logYmeas(self) -> stim.Circuit:
        "measure top and left Y logical"

        measCircuit = stim.Circuit()
        yMeas = []
        print("using topleft Ymeas")
        
        measCircuit.append("MY", 0)
        measCircuit.append("MZ", [4*i for i in range(1,3)])
        measCircuit.append("MX", [4*self.dx*i for i in range(1,3)])
        measCircuit.append("TICK")

        measlen = 5
        yMeas = [ stim.target_rec(-measlen + k)
                 for k in range(measlen)] 

        measCircuit.append("OBSERVABLE_INCLUDE", yMeas, 0) 
        return measCircuit
        
    
    def reg_3_decode(self, measlog:bool) -> stim.Circuit:
        
        dec_circ = stim.Circuit()

        cx_list = [[12,11],[10,8],[9,7],[3,1],[6,4],
            [10,11],[5,8],[9,6],[3,0],[1,4],
            [12,6],[5,1],
            [11,1],[12,7],[5,6],
            [6,11],[5,0],[2,1]]
        
        spacers = [4,9,11,14,17]

        h_list = [10, 12, 5, 2, 9, 3]

        tickrec = 0
        for cidx, i in enumerate(cx_list):
            dec_circ.append("CX",[coord_to_stim(json_to_coord[str(i[0])],df=self.dx), coord_to_stim(json_to_coord[str(i[1])],df=self.dx)] )
            if cidx == spacers[tickrec]:
                dec_circ.append("TICK")
                tickrec+=1
                
        for i in h_list:
            dec_circ.append("H",[coord_to_stim(json_to_coord[str(i)],df=self.dx)])

        dec_circ.append("TICK")
        for i in chain(range(0,6),range(7,13)):
            dec_circ.append("M",[coord_to_stim(json_to_coord[str(i)],df=self.dx)])

        for i in range(12):
            dec_circ.append("DETECTOR", stim.target_rec(-1-i))
            
        if measlog:
            print("using updated logical")
            dec_circ.append("MY", [coord_to_stim(json_to_coord[str(6)],df=self.dx)])
            dec_circ.append("OBSERVABLE_INCLUDE", arg=0, targets=[stim.target_rec(-1),
                                                                  stim.target_rec(-4),  #x stabs
                                                                  stim.target_rec(-8),
                                                                  stim.target_rec(-12),  #z stabs
                                                                  stim.target_rec(-3), 
                                                                  stim.target_rec(-13),
                                                                    ])


        return dec_circ




    
    def add_errors_to_data_qubits(self, err:int) -> stim.Circuit:
        "adds errors to data qubits in the circuit with input parameters"

        errCirc = stim.Circuit()
        indices =[]
        for i in range(13):
            # convert to stim coordinates
            indices.append( coord_to_stim(json_to_coord[str(i)],df=self.dx) )

        for e_idx in range(13):
            e = (err // (4**e_idx)) % 4
            if e % 2 == 1 :
                errCirc.append(stim.Circuit(f"X_ERROR(1) {indices[e_idx]}"))
                #print("Adding X error to qubit",  indices[e_idx])
            if e // 2 == 1 :
                errCirc.append(stim.Circuit(f"Z_ERROR(1) {indices[e_idx]}"))
                #print("Adding Z error to qubit", indices[e_idx])

        errCirc.append("TICK")
        return errCirc
    

    def grow_3u5r(self) -> stim.Circuit: #finish this!!
        mainqubs = [ row * 4 * self.dx + 4 * col 
                    for row in range(3) for col in range(3)]
        
        altqubs = [4 * self.dx * (row) + 2*self.dx + 2 + 4 * col 
                   for row in range(2) for col in range(2)]
        
        self.mgrid_qubs = np.array(mainqubs)
        self.agrid_qubs = np.array(altqubs)

        unitary_encoder = stim.Circuit()
        
        #reset all d5 new qubits
        for row in range(3):
            firstanc = 4*self.dx*row + 2
            unitary_encoder.append("R", [firstanc, firstanc+4])
            if row!= 2:
                unitary_encoder.append("R", [firstanc - 2 + 2*self.dx, firstanc+2 + 2*self.dx, firstanc+6+ 2*self.dx])
        unitary_encoder.append("TICK")
        
        for row in range(3):
            firstanc = 4*self.dx*row + 2
            unitary_encoder.append("H", [firstanc, firstanc+4])
        unitary_encoder.append("TICK")

        
        unitary_encoder.append("CNOT", flatten([[q, self.qub_below(q)] for q in self.mgrid_qubs[:6] ]))
        unitary_encoder.append("CNOT", flatten([[self.qub_below(q), q] for q in self.agrid_qubs ]))
        unitary_encoder.append("I", self.mgrid_qubs[6] ) #keep this qubit for inclusion in idling
        unitary_encoder.append("TICK")

        unitary_encoder.append("CNOT", flatten([[self.qub_toleft(q), q] for q in self.mgrid_qubs[[1,2,4,5,7,8]] ]))
        unitary_encoder.append("CNOT", flatten([[q, self.qub_toright(q)] for q in self.agrid_qubs ]))
        unitary_encoder.append("TICK")

        return unitary_encoder
    
    def __post_init__(self):
        "initialize arrays of reg3 "

        multf = 2 - int(self.smallsc)
        mainqubs = [ row * 2 * multf * self.dx + 2 * multf * col 
                    for row in range(3) for col in range(3)]
        
        if not self.smallsc:
            altqubs = [4 * self.dx * (row) + 2*self.dx + 2 + 4 * col 
                   for row in range(2) for col in range(2)]
        else:
            altqubs = [2 * (self.dx+1) * (row) + 2*(self.dx+1) + 3 + 2 * col 
                   for row in range(2) for col in range(2)]
            
        if not self.smallsc:
            bdry_as = np.array([6, 2 * self.dx, 6 * self.dx + 8, 8 * self.dx + 2])
        else:
            bdry_as = np.array([5, 2 * (self.dx +1) + 1, 4 * (self.dx +1) + 7, 6 * (self.dx +1) + 3])
        
        self.mgrid_qubs = np.array(mainqubs)
        self.agrid_qubs = np.array(altqubs)
        self.bdry_ancillas = np.array(bdry_as)


    
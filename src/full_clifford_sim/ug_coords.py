import itertools
from typing import Optional

def add_tuple(x: tuple, y: Optional[tuple] = None):
    if y is None:
        return x
    return (x[0] + y[0], x[1] + y[1])

def subtract_tuple(x: tuple, y: Optional[tuple] = None):
    if y is None:
        return x
    return (x[0] - y[0], x[1] - y[1])

def offset_loc_list(lst: list, offset: Optional[tuple] = None):
    return [add_tuple(i, offset) for i in lst]

def scale(x: tuple, f: float):
    return tuple([f * _x for _x in x])

# Rot code qubits
def get_rot_data_locs(d, a, offset: tuple = None):
    return offset_loc_list([(i*a, j*a) for j in range(d) for i in range(d)], offset)

def get_rot_init_in_x_locs(d, a, offset: tuple = None):
     ans = [(i*a, j*a) for i in range(d) for j in range(i, d)]
     assert ans[0] == (0, 0), f'ans[0] should be (0, 0) before offest, but not {ans[0]}.'
     ans = ans[1:]
     return offset_loc_list(ans, offset)

def get_rot_init_in_z_locs(d, a, offset: tuple = None):
    return offset_loc_list([(i*a, j*a) for i in range(1, d) for j in range(0, i)], offset)

def get_rot_x_anc_locs(d, a, offset: tuple = None):
    return offset_loc_list([(i*a + a/2, j*a - a/2) for i in range(d-1) for j in range(d+1) if (i+j)%2 == 0], offset)

def get_rot_z_anc_locs(d, a, offset: tuple = None):
    return offset_loc_list([(i*a - a/2, j*a + a/2) for i in range(d+1) for j in range(d-1) if (i+j)%2 == 1], offset)

def get_rot_boundary_x_anc_locs(d, a, offset: tuple = None):
    return offset_loc_list([(i*a + a/2, j*a - a/2) for i in range(d-1) for j in [0, d] if (i+j)%2 == 0], offset)

def get_rot_boundary_z_anc_locs(d, a, offset: tuple = None):
    return offset_loc_list([(i*a - a/2, j*a + a/2) for i in [0, d] for j in range(d-1) if (i+j)%2 == 1], offset)

def get_rot_bulk_x_anc_locs(d, a, offset: tuple = None):
    return offset_loc_list([(i*a + a/2, j*a - a/2) for i in range(d-1) for j in range(1, d) if (i+j)%2 == 0], offset)

def get_rot_bulk_z_anc_locs(d, a, offset: tuple = None):
    return offset_loc_list([(i*a - a/2, j*a + a/2) for i in range(1, d) for j in range(d-1) if (i+j)%2 == 1], offset)

def get_rot_bulk_anc_locs(d, a, offset: tuple = None):
    return get_rot_bulk_x_anc_locs(d, a, offset) + get_rot_bulk_z_anc_locs(d, a, offset)

def get_rot_anc_locs(d, a, offset: tuple = None):
    return get_rot_x_anc_locs(d, a, offset) + get_rot_z_anc_locs(d, a, offset)

def get_rot_boundary_anc_locs(d, a, offset: tuple = None):
    return  get_rot_boundary_x_anc_locs(d, a, offset) + get_rot_boundary_z_anc_locs(d, a, offset)

def get_rot_all_locs(d, a, offset: tuple = None):
    return get_rot_data_locs(d, a, offset) + get_rot_anc_locs(d, a, offset)

# Reg code qubits
def get_reg_data_locs(d, a, offset: tuple = None):
    lst = [(i*a, j*a) for i in range(d) for j in range(d)]
    lst += [(i*a + a/2, j*a + a/2) for i in range(d-1) for j in range(d-1)]
    return offset_loc_list(lst, offset)

def get_reg_x_anc_locs(d, a, offset: tuple = None):
    return offset_loc_list([(i*a + a/2, j*a) for i in range(d-1) for j in range(d)], offset)

def get_reg_z_anc_locs(d, a, offset: tuple = None):
    return offset_loc_list([(i*a, j*a + a/2) for i in range(d) for j in range(d-1)], offset)

def get_reg_anc_locs(d, a, offset: tuple = None):
    return get_reg_x_anc_locs(d, a, offset) + get_reg_z_anc_locs(d, a, offset)

# Gate seqs
def get_rot_gate_seq(d: int, 
                     a: int, 
                     loc_to_idx: Optional[dict] = None, 
                     flattened: bool = True, 
                     bulk_only: bool = False,
                     offset: tuple = None):
    data_locs = get_rot_data_locs(d, a, offset)
    x_anc_locs = get_rot_bulk_x_anc_locs(d, a, offset) if bulk_only else get_rot_x_anc_locs(d, a, offset)
    z_anc_locs = get_rot_bulk_z_anc_locs(d, a, offset) if bulk_only else get_rot_z_anc_locs(d, a, offset)
    anc_locs = get_rot_anc_locs(d, a, offset)
    
    # These are in unit of a/2
    x_gate_dirs = [(-1, 1), ( 1,  1), (-1, -1), (1, -1)]
    z_gate_dirs = [(-1, 1), (-1, -1), ( 1,  1), (1, -1)]

    gate_seqs = [[] for _ in range(4)]
    for i in range(4):
        for loc in x_anc_locs:
            target = add_tuple(loc, scale(x_gate_dirs[i], a/2))
            if not target in data_locs: 
                continue
            if loc_to_idx:
                loc, target = loc_to_idx[loc], loc_to_idx[target]
            gate_seqs[i].append((loc, target))
        for loc in z_anc_locs:
            control = add_tuple(loc, scale(z_gate_dirs[i], a/2))
            if not control in data_locs:
                continue
            if loc_to_idx:
                control, loc = loc_to_idx[control], loc_to_idx[loc]
            gate_seqs[i].append((control, loc))

    if not flattened:
        return gate_seqs
    gate_seqs = [list(itertools.chain(*seq)) for seq in gate_seqs]
    return gate_seqs

def get_reg_gate_seq(d: int, 
                     a: int, 
                     loc_to_idx: Optional[dict] = None, 
                     flattened: bool = True, 
                     offset: tuple = None):
    data_locs = get_reg_data_locs(d, a, offset)
    x_anc_locs = get_reg_x_anc_locs(d, a, offset)
    z_anc_locs = get_reg_z_anc_locs(d, a, offset)
    
    # These are in unit of a/2
    x_gate_dirs = [(0, 1), ( 1, 0), (-1, 0), (0, -1)]
    z_gate_dirs = [(0, 1), (-1, 0), ( 1, 0), (0, -1)]

    gate_seqs = [[] for _ in range(4)]
    for i in range(4):
        for loc in x_anc_locs:
            target = add_tuple(loc, scale(x_gate_dirs[i], a/2))
            if not target in data_locs: 
                continue
            if loc_to_idx:
                loc, target = loc_to_idx[loc], loc_to_idx[target]
            gate_seqs[i].append((loc, target))
        for loc in z_anc_locs:
            control = add_tuple(loc, scale(z_gate_dirs[i], a/2))
            if not control in data_locs:
                continue
            if loc_to_idx:
                control, loc = loc_to_idx[control], loc_to_idx[loc]
            gate_seqs[i].append((control, loc))
    return [list(itertools.chain(*seq)) for seq in gate_seqs] if flattened else gate_seqs


# Gauge fixing operators
def get_x_gauge_locs(d: int, a: int, loc: tuple, offset: tuple = None):
    # X op connecting an Z stab @ loc to top boundary
    assert loc in get_rot_z_anc_locs(d, a, offset=offset)
    loc = subtract_tuple(loc, offset)
    d, a = int(d), int(a)
    i0, j0 = add_tuple(loc, (-a/2, a/2)) if loc[0] > 0 else add_tuple(loc, (a/2, a/2))
    assert j0.is_integer()
    j0 = int(j0)
    gauge = [(i0, j * a) for j in range(j0 // a, d)]
    return offset_loc_list(gauge, offset)

def get_z_gauge_locs(d: int, a: int, loc: tuple, offset: tuple = None):
    # Z op connecting an X stab @ loc to right boundary
    assert loc in get_rot_x_anc_locs(d, a, offset=offset)
    loc = subtract_tuple(loc, offset)
    d, a = int(d), int(a)
    i0, j0 = add_tuple(loc, (a/2, -a/2)) if loc[1] > 0 else add_tuple(loc, (a/2, a/2))
    assert i0.is_integer()
    i0 = int(i0)
    gauge = [(i * a, j0) for i in range(i0 // a, d)]
    return offset_loc_list(gauge, offset)

# Other
def remove_trivial_cxs(gates: list, init_dict: dict, untouched_data: set):
    # Note: CX only
    if not untouched_data:
        return gates
    truncated_gates = []
    for c, t in zip(gates[::2], gates[1::2]):
        if (c in untouched_data or t in untouched_data) and (c in init_dict['Z'] or t in init_dict['X']):
            continue
        truncated_gates += [c, t]
    return truncated_gates


if __name__ == '__main__':
    pass
#HbondWire.f90 Python version
import MDAnalysis as mda
from MDAnalysis.analysis.distances import distance_array
import numpy as np

def get_switching(x, power_num, power_den):
    x = np.asarray(x)
    num = 1.0 - x**power_num
    den = 1.0 - x**power_den
    valid = np.abs(den) > 1e-10
    return np.where(valid, num / np.where(valid, den, 1.0), power_num / power_den)

def compute_hbond_watwater(oxy_coords, h_coords, current_box, Natoms, oxy_indices):
    Hb = np.zeros((Natoms, Natoms))
    N_oxy = len(oxy_coords)
    
    dist_OO = distance_array(oxy_coords, oxy_coords, box=current_box)
    dist_OH = distance_array(oxy_coords, h_coords, box=current_box)

    for i_idx in range(N_oxy):
        for j_idx in range(i_idx + 1, N_oxy):
            mod_rOO = dist_OO[i_idx, j_idx]
            rOiH = dist_OH[i_idx]
            rOjH = dist_OH[j_idx]

            p_angles = get_switching((rOiH + rOjH - mod_rOO) / 0.6, 8.0, 12.0)
            
            mask = rOiH < rOjH
            sum_ij = np.sum(p_angles[mask])
            sum_ji = np.sum(p_angles[~mask])

            p_oo = get_switching((mod_rOO - 2.7) / 0.5, 10.0, 20.0)

            Hb[oxy_indices[i_idx], oxy_indices[j_idx]] = sum_ij * p_oo
            Hb[oxy_indices[j_idx], oxy_indices[i_idx]] = sum_ji * p_oo

    return Hb

def compute_hbond_endpoint(p_coord, oxy_coords, h_coords, current_box, Natoms, oxy_indices, oo_threshold, own_H_coords=None):
    Hb_ep = np.zeros(Natoms)
    Hb_ep_rev = np.zeros(Natoms)
    
    dist_OO = distance_array(p_coord, oxy_coords, box=current_box)[0]
    dist_OH_p = distance_array(p_coord, h_coords, box=current_box)[0]
    dist_OH_wat = distance_array(oxy_coords, h_coords, box=current_box)

    for j_idx in range(len(oxy_coords)):
        mod_rOO = dist_OO[j_idx]
        rOiH = dist_OH_p
        rOjH = dist_OH_wat[j_idx]

        dist_mask = (rOiH < 1.3) | (rOjH < 1.3)
        rOiH_f = rOiH[dist_mask]
        rOjH_f = rOjH[dist_mask]

        p_angles = get_switching((rOiH_f + rOjH_f - mod_rOO) / 0.6, 8.0, 12.0)
        mask = rOiH_f < rOjH_f
        sum_pj = np.sum(p_angles[mask])
        sum_jp = np.sum(p_angles[~mask])

        if own_H_coords is not None and len(own_H_coords) > 0:
            dist_pH_own = distance_array(p_coord, own_H_coords, box=current_box)[0]
            dist_jH_own = distance_array(oxy_coords[j_idx:j_idx+1], own_H_coords, box=current_box)[0]

            dist_mask_own = (dist_pH_own < 1.3) | (dist_jH_own < 1.3)
            rOiH_own = dist_pH_own[dist_mask_own]
            rOjH_own = dist_jH_own[dist_mask_own]

            p_angles_own = get_switching((rOiH_own + rOjH_own - mod_rOO) / 0.6, 8.0, 12.0)
            mask_own = rOiH_own < rOjH_own
            sum_pj += np.sum(p_angles_own[mask_own])
            sum_jp += np.sum(p_angles_own[~mask_own])

        p_oo = get_switching((mod_rOO - oo_threshold) / 0.5, 10.0, 20.0)

        Hb_ep[oxy_indices[j_idx]] = sum_pj * p_oo
        Hb_ep_rev[oxy_indices[j_idx]] = sum_jp * p_oo

    return Hb_ep, Hb_ep_rev

def analyse(topfile, trajfile, Nsolute, point, box):
    u = mda.Universe(topfile, trajfile)
    p1, p2 = point[0] - 1, point[1] - 1
    Nwat = Nsolute + 1
    atom = np.array([a.name.strip() for a in u.atoms])
    Natoms = len(u.atoms)

    oxy_indices = list(range(Nwat - 1, Natoms, 3))
    h_indices = [k for k in range(Nwat - 1, Natoms) if atom[k] == "H"]
    own_H1_indices = [k for k in range(p1 + 1, p1 + 3) if atom[k] == "H"]

    output_files = {
        "nuc": open("nuc_wat_chains", "w"),
        "chain1": open("1_water_chains", "w"),
        "chain2": open("2_water_chains", "w"),
        "chain3": open("3_water_chains", "w"),
    }
    THRESHOLD = 0.8

    for f_idx, ts in enumerate(u.trajectory, start=1):
        X = u.atoms.positions
        current_box = ts.dimensions if box is None else np.array(box)

        oxy_coords = X[oxy_indices]
        h_coords = X[h_indices]
        p1_coord = X[p1:p1+1]
        p2_coord = X[p2:p2+1]
        own_H1_coords = X[own_H1_indices] if own_H1_indices else None

        Hb = compute_hbond_watwater(oxy_coords, h_coords, current_box, Natoms, oxy_indices)

        Hb_p1, Hb_p1_rev = compute_hbond_endpoint(
            p1_coord, oxy_coords, h_coords, current_box, Natoms, oxy_indices, 
            oo_threshold=2.7, own_H_coords=own_H1_coords
        )

        Hb_p2, Hb_p2_rev = compute_hbond_endpoint(
            p2_coord, oxy_coords, h_coords, current_box, Natoms, oxy_indices, 
            oo_threshold=3.7
        )

        for j in oxy_indices:
            Hb[p1, j], Hb[j, p1] = Hb_p1[j], Hb_p1_rev[j]
            Hb[p2, j], Hb[j, p2] = Hb_p2[j], Hb_p2_rev[j]

        if Hb[p1, p2] >= THRESHOLD:
            output_files["nuc"].write(f"{f_idx} {p2+1} {Hb[p1,p2]:.6f}\n")

        for i in oxy_indices:
            if Hb[p1, i] * Hb[i, p2] >= THRESHOLD:
                output_files["chain1"].write(f"{f_idx} {i+1} {Hb[p1,i]:.6f} {Hb[i,p2]:.6f}\n")
            for j in oxy_indices:
                if Hb[p1, i] * Hb[i, j] * Hb[j, p2] >= THRESHOLD:
                    output_files["chain2"].write(f"{f_idx} {i+1} {j+1} {Hb[p1,i]:.6f} {Hb[i,j]:.6f} {Hb[j,p2]:.6f}\n")
                for k in oxy_indices:
                    if Hb[p1, i] * Hb[i, j] * Hb[j, k] * Hb[k, p2] >= THRESHOLD:
                        output_files["chain3"].write(f"{f_idx} {i+1} {j+1} {k+1} {Hb[p1,i]:.6f} {Hb[i,j]:.6f} {Hb[j,k]:.6f} {Hb[k,p2]:.6f}\n")

        if f_idx % 100 == 0:
            print(f"Processed frame {f_idx}")

    print(f"Last frame: {f_idx}")
    for fh in output_files.values():
        fh.close()

if __name__ == "__main__":
    analyse(
        topfile="topology.prmtop",
        trajfile="trajectory.dcd",
        Nsolute=3,               
        point=[1, 4],          
        box=None            
    )

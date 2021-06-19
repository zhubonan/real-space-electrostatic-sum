
#cimport scipy.special.cython_special as cython_special
"""
Main module for summation
Based on real_space_eletrostatic_sum.f90
"""
from itertools import product
from math import pi, sqrt, erfc, erf, exp
import numpy as np

cpi = pi
sqrt_pi = sqrt(pi)
one_third = 1.0 / 3.0

def nacl_base():
    """Base case for Al test"""
    # lattice vectors
    a_1 = np.array([1.0, 1.0, 0.0])
    a_2 = np.array([0.0, 1.0, 1.0])
    a_3 = np.array([1.0, 0.0, 1.0])
    lattice = np.stack([a_1, a_2, a_3], axis=0)

    # length scale and cutoff
    h_max = np.sqrt(4.0 / 3.0)
    r_d_hat = 3.0

    # ionic positions and charge
    positions = np.zeros([2, 3], dtype=np.float)
    positions[1, :] = [1., 1., 1.]
    chg = np.array([-1., 1.])

    rc = 3.0 * r_d_hat * 2 * h_max
    rd = r_d_hat * h_max

    ewald = -1.747564594633
    return lattice, positions, chg, rc, rd, ewald

nacl_param = nacl_base()

def energy(lattice, positions, z, rc, rd):
    """
    Compute the energies

    Args:
      lattice: a (3,3) array of row vectors
      postions: a (N, 3) array of the positions
      z: a (N) array of the charges 
      rc: cut off radius
      rd: adaptive cut off

    Retruns:
      A tuple of energy, E_i matrix, delta_E_i array
    """

    a1, a2, a3 = lattice
    vol = abs(np.cross(a1, a2) @ a3)
    rho_positive = z[z > 0].sum() / vol  # Neutralising backgroud for positive 
    rho_negative = z[z < 0].sum() / vol  # Neutralising backgroud for negative 
    nions = positions.shape[0]
    assert len(z) == nions, "Mismatch in the ionic positions and the charges array"

    # Obtain the reciprocal lattice vectors
    invl = np.linalg.inv(lattice)

    # Distances between the planes 
    # d_100, d_010, d_001 = 1.0 / np.linalg.norm(invl, axis=1)

    # Compute the maximum shifts needed
    shift1max, shift2max, shift3max = np.ceil(rc * np.linalg.norm(invl, axis=1)).astype(int)

    # Loop over the cells and pre-compute the shift vectors
    nshifts = (shift1max * 2 + 1) * (shift2max * 2 + 1) * (shift3max * 2 + 1)

    # Pre-allocate the shift vectors
    shift_vectors = np.zeros((nshifts, 3), dtype=np.float)
    shift_indices = np.zeros((nshifts, 3), dtype=np.int)
    iter_tmp = product(range(-shift3max, shift3max+1),
                                            range(-shift2max, shift2max+1),
                                            range(-shift1max, shift1max+1),)
    itmp = 0
    for shift3, shift2, shift1 in iter_tmp:
        shift_vectors[itmp, :] = shift1 * a1 + shift2 * a2 + shift3 * a3
        shift_indices[itmp, :] = (shift1, shift2, shift3)
        itmp += 1

    # Pre-allocate E_i matrix and \Delta E_i matrix
    E_i = np.zeros((nions, nions), dtype=np.float)
    delta_E_i = np.zeros(nions, dtype=np.float)


    energy = 0.0

    for i in range(nions):

        # Prepare for loop over neighboring ions
        ei = 0.0
        qi_positive = z[i]
        qi_negative = z[i]

        # Loop over the cells
        iter_tmp = product(range(-shift3max, shift3max+1),
                                              range(-shift2max, shift2max+1),
                                              range(-shift1max, shift1max+1),)
        for origin_j, shift_idx in zip(shift_vectors, shift_indices):
           
            # Loop over the other ions
            for j in range(nions):
               
                # Connect interactive with itself
                if (i==j) and all(shift_idx == 0):
                    continue

                # Distance
                rij = comp_distance(positions[i], positions[j], origin_j)

                # No need to continue if it is outside the cut-off
                if rij > rc:
                    continue

                # Update energy
                E_i[i, j] = erfc(rij / rd) / rij

                # Update the total change inside the cut-off sphere (Q_i)
                # Positive and negative terms have to be counted separately
                if z[j] > 0.:
                    qi_positive += z[j]
                if z[j] < 0.:
                    qi_negative += z[j]

        # Apply 1/2 z[i] factor to energy
        E_i[i, :] = E_i[i, :] * 0.5 * z[i]
        ei = E_i[i, :].sum()

        # Compute the correction term for positive and negative ions separately
        # Correction terms - positive
        if rho_positive != 0.0:
            ra = (3.0 * qi_positive / (4.0 * pi * rho_positive)) ** one_third  # adaptive cutoff
            delta_ei_positive = comp_delta_ei(ra, rho_positive, z[i], rd)
        else:
            delta_ei_positive = 0.0

        # Correction terms - negative
        if rho_negative != 0.0:
            ra = (3.0 * qi_negative / (4.0 * pi * rho_negative)) ** one_third  # adaptive cutoff
            # Accumulate the correction term
            delta_ei_negative = comp_delta_ei(ra, rho_negative, z[i], rd)
        else:
            delta_ei_negative = 0.0


        # Compute the total correction term
        delta_ei = delta_ei_negative + delta_ei_positive
        delta_E_i[i] = delta_ei
        energy += ei + delta_ei

    return float(energy), E_i, delta_E_i


def comp_delta_ei(ra, rho, zi, rd):
    """
    Compute the correction term - Eq(19) in the reference
    """
    value = - pi * zi * rho * ra * ra + pi * zi * rho * ( ra * ra - rd * rd / 2.0)\
                    * erf(ra / rd) + sqrt_pi * zi * rho * ra * rd * exp(-ra * ra / (rd * rd)) \
                    - 1.0 / (sqrt_pi * rd) * zi * zi
    return value

def comp_distance(v1, v2, shift):
    """
    Compute the distance between v1 and v2 subject to a shift vector
    """
    vtmp = v1 - (v2 + shift)
    return np.linalg.norm(vtmp)


energy(*nacl_param[:-1])

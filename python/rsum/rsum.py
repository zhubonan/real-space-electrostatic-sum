"""
Main module for summation
Based on real_space_eletrostatic_sum.f90
"""
from itertools import product
from math import pi, sqrt, erfc, erf, exp
import numpy as np

sqrt_pi = sqrt(pi)
one_third = 1.0 / 3.0

def energy(lattice, positions, z, rc, rd):
    """
    Compute the energies

    Args:
      lattice: a (3,3) array of row vectors
      postions: a (N, 3) array of the positions
      z: a (N) array of the charges 
      rc: cut off radius
      rd: adaptive cut off
    """

    a1, a2, a3 = lattice
    vol = abs(np.cross(a1, a2) @ a3)
    rho = z.sum() / vol   # Average charge density
    nions = positions.shape[0]
    assert len(z) == nions, "Mismatch in the ionic positions and the charges array"

    # Obtain the reciprocal lattice vectors
    invl = np.linalg.inv(lattice)

    # Distances between the planes 
    # d_100, d_010, d_001 = 1.0 / np.linalg.norm(invl, axis=1)

    # Compute the maximum shifts needed
    shift1max, shift2max, shift3max = np.ceil(rc * np.linalg.norm(invl, axis=1)).astype(int)


    energy = 0.0

    for i in range(nions):

        # Prepare for loop over neighboring ions
        ei = 0.0
        qi = z[i]

        # Loop over the cells
        iter_tmp = product(range(-shift3max, shift3max+1),
                                              range(-shift2max, shift2max+1),
                                              range(-shift1max, shift1max+1),)
        for shift3, shift2, shift1 in iter_tmp:
           
            origin_j = shift1 * a1 + shift2 * a2 + shift3 * a3

            # Loop over the other ions
            for j in range(nions):
               
                # Connect interactive with itself
                if (i==j) and all([shift1==0, shift2==0, shift3==0]):
                    continue

                # Distance
                vij = positions[i] - (positions[j] + origin_j)
                rij = sqrt(np.dot(vij, vij))

                # No need to continue if it is outside the cut-off
                if rij > rc:
                    continue

                # Update energy and charge
                ei = ei + z[j] * erfc(rij / rd) / rij
                qi = qi + z[j]

        # Apply 1/2 z[i] factor to energy
        ei = 0.5 * z[i] * ei

        # Correction terms
        ra = (3.0 * qi / (4.0 * pi * rho)) ** one_third  # adaptive cutoff
        ei = ei - pi * z[i] * rho * ra * ra + pi * z[i] * rho * ( ra * ra - rd * rd / 2.0)\
            * erf(ra / rd) + sqrt_pi * z[i] * rho * ra * rd * exp(-ra * ra / (rd * rd)) \
            - 1.0 / (sqrt_pi * rd) * z[i] * z[i]

        energy += ei

    return float(energy)
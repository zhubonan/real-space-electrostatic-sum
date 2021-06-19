"""
Test for the rsum module
"""

import pytest
import numpy as np
from rsum.rsum import energy


def _al_base():
    """Base case for Al test"""
    a_1 = np.array([5.41141973394663, 0.00000000000000, 0.00000000000000])
    a_2 = np.array([2.70570986697332, 4.68642696013821, 0.00000000000000])
    a_3 = np.array([2.70570986697332, 1.56214232004608, 4.41840571073226])
    lattice = np.stack([a_1, a_2, a_3], axis=0)
    positions = np.zeros([1, 3])
    chg = 3.0 * np.ones(positions.shape[0])
    h_max = 4.42
    r_d_hat = 2.0

    rc = 3.0 * r_d_hat * 2 * h_max
    rd = r_d_hat * h_max

    ewald = -2.69595457432924945
    return lattice, positions, chg, rc, rd, ewald


def _nacl_base():
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


al_base = pytest.fixture(_al_base)
nacl_base = pytest.fixture(_nacl_base)


def test_energy_al(al_base):
    """
    Test the computeation of energy
    """

    eng, E_i, delta_E_i = energy(*al_base[:5])
    ewald = al_base[-1]
    assert eng == E_i[0, 0] + delta_E_i[0]
    assert pytest.approx(eng == ewald, rel=1e-9)


def test_energy_nacl(nacl_base):
    """
    Test the computeation of energy
    """

    eng, E_i, delta_E_i = energy(*nacl_base[:5])
    ewald = nacl_base[-1]
    assert eng == E_i.sum() + delta_E_i.sum()
    assert pytest.approx(eng == ewald, rel=1e-9)
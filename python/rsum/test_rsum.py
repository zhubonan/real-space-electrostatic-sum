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

    rc = 3.0 * r_d_hat ** 2 * h_max
    rd = r_d_hat * h_max

    ewald = -2.69595457432924945
    return lattice, positions, chg, rc, rd, ewald

def _si_base():
    """Base case for Al test"""
    # Si
    a_1 = np.array([7.25654832321381, 0.00000000000000, 0.00000000000000])
    a_2 = np.array([3.62827416160690, 6.28435519169252, 0.00000000000000])
    a_3 = np.array([3.62827416160690, 2.09478506389751, 5.92494689524090])
    lattice = np.stack([a_1, a_2, a_3], axis=0)
    loc = np.array([[0.0,  0.0,  0.0],
                    [0.25, 0.25, 0.25]])
    positions = (np.vstack((a_1, a_2, a_3)).T).dot(loc.T).T # to cartesian
    chg = 4.0 * np.ones(positions.shape[0])
    h_max = 5.92
    r_d_hat = 2.0
    ewald = -8.39857465282205418
    rc = 3.0 * r_d_hat ** 2 * h_max
    rd = r_d_hat * h_max

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

    rc = 3.0 * r_d_hat ** 2 * h_max
    rd = r_d_hat * h_max

    ewald = -1.747564594633
    return lattice, positions, chg, rc, rd, ewald


al_base = pytest.fixture(_al_base)
nacl_base = pytest.fixture(_nacl_base)
si_base = pytest.fixture(_si_base)


def test_energy_al(al_base):
    """
    Test the computeation of energy
    """

    eng, E_i, delta_E_i = energy(*al_base[:5])
    ewald = al_base[-1]
    assert eng == E_i[0, 0] + delta_E_i[0]
    assert eng == pytest.approx(ewald, abs=1e-9)


def test_energy_si(si_base):
    """
    Test the computeation of energy
    """

    eng, E_i, delta_E_i = energy(*si_base[:5])
    ewald = si_base[-1]
    assert eng == pytest.approx(ewald, abs=1e-9)


def test_energy_nacl(nacl_base):
    """
    Test the computeation of energy
    """

    eng, E_i, delta_E_i = energy(*nacl_base[:5])
    ewald = nacl_base[-1]
    assert eng == E_i.sum() + delta_E_i.sum()
    assert eng == pytest.approx(ewald, rel=1e-9)
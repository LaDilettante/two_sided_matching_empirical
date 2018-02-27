import numpy as np
from numpy.random import RandomState
import pytest

import abm

rng1 = RandomState(1)
@pytest.mark.parametrize("er1, er2, ee", [
    (abm.Employer(w=np.array([100, 100]), beta=np.array([-99]), rng=rng1),
     abm.Employer(w=np.array([1, 1]), beta=np.array([-99]), rng=rng1),
     abm.Employee(x=np.array([-99]), alpha=np.array([0.5, 0.6]), rng=rng1)),
])
def test_gumbel_rvs_works(ee, er1, er2):
    ee.pick_best_employer([er1, er2], np.array([0, 1])) == 1
    assert np.allclose(ee.wa, np.array([110.13397001, 2.21457862]))

@pytest.mark.parametrize("er, ee", [
    (abm.Employer(w=np.array([-99]), beta=np.array([-10, 1, 1])),
     abm.Employee(x=np.array([1, 2, 2]), alpha=np.array([-99]))),
])
def test_employer_make_offer(er, ee):
    assert er.make_offer([ee]) == [0] 


@pytest.mark.parametrize("er1, er2, ee", [
    (abm.Employer(w=np.array([100, 100]), beta=np.array([-99])),
     abm.Employer(w=np.array([1, 1]), beta=np.array([-99])),
     abm.Employee(x=np.array([-99]), alpha=np.array([0.5, 0.5]))),
])
def test_employee_pick_offer(ee, er1, er2):
    assert ee.pick_best_employer([er1, er2], np.array([0, 1])) == 1

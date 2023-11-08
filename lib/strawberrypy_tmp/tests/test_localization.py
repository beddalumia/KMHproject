import numpy as np

import os, sys
currentdir = os.path.dirname(os.path.realpath(__file__))
parentdir = os.path.dirname(currentdir)
sys.path.append(parentdir)

from strawberrypy_tmp import FiniteModel
from strawberrypy_tmp.example_models import kane_mele_tbmodels, kane_mele_pythtb

def test_localization_marker():
    # Construction of the Kane-Mele model with TBmodels and PythTB
    kmmodel_tbm = kane_mele_tbmodels(rashba = 1, esite = 3, spin_orb = 0.3, L = 1)
    kmmodel_pythtb = kane_mele_pythtb(rashba = 1, esite = 3, spin_orb = 0.3, L = 1)

    # Construct the two finite models
    kmmodel_tbm = FiniteModel(tbmodel = kmmodel_tbm, nx_sites = 10, ny_sites = 10, spinful = True, mode = 'tb')
    kmmodel_pythtb = FiniteModel(tbmodel = kmmodel_pythtb, nx_sites = 10, ny_sites = 10, spinful = True, mode = 'tb')

    # Compute the local Chern marker on the lattice
    loc_tbm = kmmodel_tbm.localization_marker()
    loc_pythtb = kmmodel_pythtb.localization_marker()

    # Check the two models give the same results
    assert np.allclose(loc_tbm, loc_pythtb)

    # Add Anderson disorder
    seed = 938418523
    kmmodel_tbm.add_onsite_disorder(w = 0, seed = seed)
    kmmodel_pythtb.add_onsite_disorder(w = 0, seed = seed)

    # Add random vacancies to the system by removing a lattice site
    random_vacancies = [[np.random.randint(10), np.random.randint(10), np.random.randint(2)] for _ in range(5)]
    kmmodel_tbm.add_vacancies(vacancies_list = random_vacancies)
    kmmodel_pythtb.add_vacancies(vacancies_list = random_vacancies)

    # Evaluate the local Chern marker averaging over a small region of the sample
    loc_tbm = kmmodel_tbm.localization_marker(direction = 0, start = 5, macroscopic_average = True, cutoff = 2)
    loc_pythtb = kmmodel_pythtb.localization_marker(direction = 0, start = 5, macroscopic_average = True, cutoff = 2)

    # Check the two models give the same results
    assert np.allclose(loc_tbm, loc_pythtb)
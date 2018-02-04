import numpy as np
from itertools import compress
import rpy2.robjects as ro
from rpy2.robjects import numpy2ri
numpy2ri.activate()
from rpy2.robjects.packages import importr

class Employer:
    '''
    w : p_j vector of characteristics
    beta : p_i vector of preferences for employee's characteristics
    '''
    def __init__(self, w=None, beta=None, reserve_utility=0):
        self.w = np.random.normal(size=2) if w is None else w
        self.beta = np.random.normal(size=3) if beta is None else w
        self.reserve_utility = reserve_utility

    def make_offer(self, employee_list):

        return [1 if self.beta.dot(employee.x) > self.reserve_utility else 0
                for employee in employee_list]

class Employee:
    '''
    x : p_i vector of characteristics
    alpha : p_j vector of preferences for employer's characteristics
    '''

    alpha = np.random.normal(size=2)

    def __init__(self, x=None, alpha=None):
        self.x = np.random.normal(size=3) if x is None else x

    def offer_list(self, offer_logical, employer_list):
        return employer_list

    def pick_best_employer(self, offer_list):
        employer_values = np.array([self.alpha.dot(employer.w)
                                    for employer in offer_list])
        return np.argmax(employer_values)

class Model:

    def __init__(self, n_employer, n_employee):
        self.num_employer = n_employer
        self.num_employee = n_employee

    def matching_process(self):
        employer_list = [Employer() for i in range(self.num_employer)]
        employee_list = [Employee() for i in range(self.num_employee)]

        # Employer making offers
        opportunity_set = np.transpose(np.array([employer.make_offer(employee_list)
                                                for employer in employer_list]))
        opportunity_set[:, 0] = 1 # Unemployement is always available

        # Employee choosing where to work
        choice = np.array([employee.pick_best_employer(compress(employer_list,
                                                                opportunity_set[idx, :]))
                           for idx, employee in enumerate(employee_list)])

        # Package the characteristics for R
        xx = np.array([e.x for e in employee_list])
        ww = np.array([e.w for e in employer_list])

        # Package the preferences for R
        alpha = np.unique(np.array([e.alpha for e in employee_list]), axis=0)
        assert alpha.shape[0] == 1 # There is only one set of alpha
        beta = np.array([e.beta for e in employer_list])

        return (alpha, beta, ww, xx, opportunity_set, choice)

my_model = Model(n_employer=5, n_employee=20)
alpha, beta, ww, xx, opp, choice = my_model.matching_process()

importr("MASS", lib_loc="/usr/lib/R/library/")
importr("lattice", lib_loc="/usr/lib/R/library/")
importr("Matrix", lib_loc="/usr/lib/R/library/")

match2sided = ro.r['match2sided']
match2sided(iter=200, eps_alpha = 0.4, eps_beta = 0.05, 
                  frac_beta = 1, frac_opp = 0.25,
            ww = ww, xx = xx, choice = choice, opp = opp)

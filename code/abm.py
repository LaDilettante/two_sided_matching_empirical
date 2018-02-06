import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from itertools import compress
import rpy2.robjects as ro
from rpy2.robjects import numpy2ri
numpy2ri.activate()
from rpy2.robjects.packages import importr
from rpy2.robjects.packages import SignatureTranslatedAnonymousPackage

class Employer:
    '''
    w : p_j vector of characteristics
    beta : p_i vector of preferences for employee's characteristics
    '''
    def __init__(self, w=None, beta=None, reserve_utility=0):
        self.w = np.random.normal(size=2) if w is None else w
        self.beta = np.random.normal(size=3) if beta is None else beta
        self.reserve_utility = reserve_utility

    def make_offer(self, employee_list):

        return [1 if self.beta.dot(employee.x) > self.reserve_utility else 0
                for employee in employee_list]

class Employee:
    '''
    x : p_i vector of characteristics
    alpha : p_j vector of preferences for employer's characteristics
    '''

    def __init__(self, x=None, alpha=None):
        self.x = np.random.normal(size=3) if x is None else x
        self.alpha = np.random.normal(size=2) if alpha is None else alpha

    def offer_list(self, offer_logical, employer_list):
        return employer_list

    def pick_best_employer(self, employer_list, offer_list):
        """
        offer_list : a 0/1 mask of which employer made an offer
        """
        # employer values is -inf for the ones that are not offered
        employer_values = np.array([float("-inf")] * len(employer_list))
        wa = np.array([self.alpha.dot(employer.w) for employer in employer_list])
        employer_values[offer_list.astype("bool")] = wa[offer_list.astype("bool")]
        return np.argmax(employer_values)

class Model:

    def __init__(self, n_employer, n_employee):
        self.num_employer = n_employer
        self.num_employee = n_employee

    def matching_process(self, employer_list=None, employee_list=None):
        if employer_list is None:
            employer_list = [Employer() for i in range(self.num_employer)]
        if employee_list is None:
            employee_list = [Employee() for i in range(self.num_employee)]

        # Employer making offers
        opportunity_set = np.transpose(np.array([employer.make_offer(employee_list)
                                                for employer in employer_list]))
        opportunity_set[:, 0] = 1 # Unemployment is always available

        # Employee choosing where to work
        choice = np.array([ee.pick_best_employer(employer_list, opportunity_set[idx, :])
                           for idx, ee in enumerate(employee_list)])

        # Compute the characteristics for R
        xx = np.array([e.x for e in employee_list])
        xx = np.column_stack((np.ones(xx.shape[0]), xx))
        ww = np.array([e.w for e in employer_list])

        # Compute the preferences for R
        alpha = np.unique(np.array([e.alpha for e in employee_list]), axis=0)
        assert alpha.shape[0] == 1 # There is only one set of alpha
        beta = np.array([e.beta for e in employer_list])

        # Compute the observed opportunity set for R
        observed_opp = np.zeros((self.num_employee, self.num_employer))
        observed_opp[:, 0] = 1 # Unemployement is always available
        observed_opp[np.arange(self.num_employee), choice] = 1 # Accepted jobs are available

        return (alpha, beta, ww, xx, choice, opportunity_set, observed_opp)









# my_model = Model(n_employer=10, n_employee=100)
# alpha, beta, ww, xx, choice, opp, observed_opp = my_model.matching_process()

# plt.imshow(opp) ; plt.show()
# plt.imshow(observed_opp) ; plt.show()

# np.save("ww.npy", ww)
# np.save("xx.npy", xx)
# np.save("opp.npy", observed_opp)
# np.save("choice.npy", choice)


# importr("MASS", lib_loc="/usr/lib/R/library/")
# importr("lattice", lib_loc="/usr/lib/R/library/")
# importr("Matrix", lib_loc="/usr/lib/R/library/")

# with open('match2sided.R', 'r') as myfile:
#     data = myfile.read()
#     matchpack = SignatureTranslatedAnonymousPackage(data, "matchpack")

# iter = 5000
# res = matchpack.match2sided(iter=iter,
#                             C_alpha = np.diag(np.repeat(0.5 ** 2, ww.shape[1])),
#                             C_beta = np.diag(np.repeat(0.5 ** 2, xx.shape[1])),
#                             frac_opp = 0.2,
#                             ww = ww, xx = xx,
#                             choice = choice, opp = observed_opp)

# alpha_mcmc = np.array(res.rx2('alpha'))
# print('True {}. Posterior mean {} +- sd {}'
#       .format(alpha,
#               np.mean(alpha_mcmc[:iter // 2, :], axis=0),
#               np.std(alpha_mcmc[:iter // 2, :], axis=0)))
# plt.plot(alpha_mcmc[:iter // 2, :]) ; plt.show()

import numpy as np
from numpy.random import RandomState
from scipy.stats import gumbel_r

class Employer:
    '''
    w : p_j vector of characteristics
    beta : p_i vector of preferences for employee's characteristics
    p_j : the number of employer characteristics
    '''
    def __init__(self, w=None, beta=None, p_j=None, p_i=None, rng=None):
        self.w = np.random.normal(size=p_j) if w is None else w
        self.beta = np.random.normal(size=p_i) if beta is None else beta
        self.rng = rng

    def make_offer(self, employee_list):
        '''Make 0 / 1 offer to a list of employee'''
        utility_hiring = np.array([self.beta.dot(ee.x) + gumbel_r().rvs(random_state=self.rng)
                                   for ee in employee_list])
        utility_not_hiring = np.array([gumbel_r().rvs(random_state=self.rng)
                                       for ee in employee_list])

        return (utility_hiring > utility_not_hiring).astype("float64")

class Employee:
    '''
    x : p_i vector of characteristics
    alpha : p_j vector of preferences for employer's characteristics
    '''

    def __init__(self, x=None, alpha=None, p_i=None, p_j=None, rng=None):
        self.x = np.random.normal(size=p_i) if x is None else x
        self.alpha = np.random.normal(size=p_j) if alpha is None else alpha
        self.rng = rng

    def pick_best_employer(self, employer_list, offer_list):
        """
        offer_list : a 0/1 mask of which employer made an offer
        """
        # employer values is -inf for the ones that are not offered
        utilities = np.array([float("-inf")] * len(employer_list))
        wa = np.array([self.alpha.dot(employer.w) + gumbel_r().rvs(random_state=self.rng)
                       for employer in employer_list])
        wa = np.squeeze(wa)
        self.wa = wa
        utilities[offer_list.astype("bool")] = wa[offer_list.astype("bool")]
        return np.argmax(utilities)

class Model:

    def __init__(self, employer_list, employee_list, rng=None):
        self.employer_list = employer_list
        self.employee_list = employee_list
        self.num_employer = len(employer_list)
        self.num_employee = len(employee_list)
        self.rng = rng
        
    def matching_process(self, employer_list=None, employee_list=None):
        if employer_list is None:
            employer_list = [Employer() for i in range(self.num_employer)]
        if employee_list is None:
            employee_list = [Employee() for i in range(self.num_employee)]

        # Employer making offers
        opportunity_set = np.transpose(np.array([employer.make_offer(employee_list)
                                                for employer in employer_list]))
        opportunity_set[:, 0] = 1 # Unemployment is always available
        opportunity_set = opportunity_set.astype("float64")

        # Employee choosing where to work
        choice = np.array([ee.pick_best_employer(employer_list, opportunity_set[idx, :])
                           for idx, ee in enumerate(employee_list)])

        # Compute the characteristics for R
        xx = np.array([e.x for e in employee_list])
        ww = np.array([e.w for e in employer_list])

        # Compute the preferences for R
        alpha = np.unique(np.array([e.alpha for e in employee_list]), axis=0)
        assert alpha.shape[0] == 1 # There is only one set of alpha
        beta = np.array([e.beta for e in employer_list])

        # Compute the observed opportunity set for R
        observed_opp = np.zeros((self.num_employee, self.num_employer))
        observed_opp[:, 0] = 1 # Unemployement is always available
        observed_opp[np.arange(self.num_employee), choice] = 1 # Accepted jobs are available

        return (alpha, beta, ww, xx, choice,
                opportunity_set, observed_opp)

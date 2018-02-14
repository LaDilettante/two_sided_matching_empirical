import numpy as np
from numpy.random import RandomState
from scipy.stats import gumbel_r

rng = RandomState(42)

class Employer:
    '''
    w : p_j vector of characteristics
    beta : p_i vector of preferences for employee's characteristics
    job : a function that returns a job. A job is a numpy array of characteristics
    '''
    def __init__(self, offer, w=None, beta=None):
        self.w = np.random.normal(size=2) if w is None else w
        self.beta = np.random.normal(size=3) if beta is None else beta
        self.offer = offer

    def create_offer(self, employee):
        return self.offer(self, employee)

    def offer_or_not(self, employee):
        '''Make 0 / 1 offer to a list of employee'''
        utility_hiring = self.beta.dot(employee.x) + gumbel_r().rvs(random_state=rng)
        utility_not_hiring = gumbel_r().rvs(random_state=rng)

        return (utility_hiring > utility_not_hiring).astype("float64")

class Employee:
    '''
    x : p_i vector of characteristics
    alpha : p_j vector of preferences for employer's characteristics
    '''

    def __init__(self, x=None, alpha=None):
        self.x = np.random.normal(size=3) if x is None else x
        self.alpha = np.random.normal(size=2) if alpha is None else alpha

    def pick_best_employer(self, employer_list, offer_list):
        """
        offer_list : a 0/1 mask of which employer made an offer
        """
        # employer values is -inf for the ones that are not offered
        utilities = np.array([float("-inf")] * len(employer_list))
        jobs = [er.create_offer(self) for er in employer_list]
        self.jobs = jobs
        job_utilities = np.array([self.alpha.dot(job) + gumbel_r().rvs(random_state=rng)
                                  for job in jobs])
        self.job_utilities = job_utilities
        utilities[offer_list.astype("bool")] = job_utilities[offer_list.astype("bool")]
        return np.argmax(utilities)

class Model:

    def __init__(self, employer_list, employee_list):
        self.employer_list = employer_list
        self.employee_list = employee_list

    def matching_process(self):
        '''Run the matching process'''
        num_employer = len(self.employer_list)
        num_employee = len(self.employee_list)

        # Employer making offers
        true_opp = np.array([[er.offer_or_not(ee) for er in self.employer_list]
                             for ee in self.employee_list])
        true_opp[:, 0] = 1 # Unemployment is always available
        true_opp = true_opp.astype("float64")

        # Employee choosing where to work
        choice = np.array([ee.pick_best_employer(self.employer_list, true_opp[idx, :])
                           for idx, ee in enumerate(self.employee_list)])

        # Compute the characteristics for R
        xx = np.array([ee.x for ee in self.employee_list])
        ww = np.array([er.w for er in self.employer_list])

        # Compute the job offers for R
        jobs = np.array([ee.jobs for ee in self.employee_list])
        job_utilities = np.array([ee.job_utilities for ee in self.employee_list])

        # Compute the preferences for R
        alpha = np.unique(np.array([ee.alpha for ee in self.employee_list]), axis=0)
        assert alpha.shape[0] == 1 # There is only one set of alpha
        beta = np.array([er.beta for er in self.employer_list])

        # Compute the observed opportunity set for R
        obs_opp = np.zeros((num_employee, num_employer))
        obs_opp[:, 0] = 1 # Unemployement is always available
        obs_opp[np.arange(num_employee), choice] = 1 # Accepted jobs are available

        return (alpha.astype("float64"), beta.astype("float64"),
                ww.astype("float64"), xx.astype("float64"),
                jobs, job_utilities,
                choice.astype("float64"), true_opp, obs_opp)

import numpy as np
from scipy import stats
import pandas as pd
import matplotlib.pyplot as plt

from numpy.random import RandomState

import abm

df_employee = pd.read_csv('../data_clean/labor_employee.csv', header=0)
df_employer = pd.read_csv('../data_clean/labor_employer_occ5.csv', header=0)

# ---- Prepare xx ----
df_xx = (df_employee
      [["educ", "age", "nonwhite"]])
df_xx.insert(0, "intercept", 1)

# ---- Prepare ww ----
df_ww = df_employer.loc[:, ['pres', 'aut']]

n_i = df_xx.shape[0]
p_i = df_xx.shape[1]
n_j = df_ww.shape[0]
p_j = df_ww.shape[1]

# Create employers
rng_beta = RandomState(1)
# true_beta = rng_beta.multivariate_normal(mean=[-15, 0.5, 0.2,  2],
#                                          cov=[[3, 0, 0, 0], 
#                                               [0, 0.15, 0, 0],
#                                               [0, 0, 0.1, 0],
#                                               [0, 0, 0, 0.5]],
#                                          size=n_j)
rng = RandomState(1)
true_beta = np.array([[0, 0, 0, 0],
                      [-24, 1.3, 0.1, 1],
                      [-22, 1, 0.2, 1],
                      [-9, 0.75, -0.05, 0],
                      [-8, 0.5, 0.02, 0],
                      [-6, 0.5, -0.01, 1]])
employer_list = [abm.Employer(w=df_ww.loc[j, :], beta=true_beta[j, :], rng=rng) for j in range(n_j)]

# Create employees
true_alpha = np.array([0.1, 1])
employee_list = [abm.Employee(x=df_xx.loc[i, :], alpha=true_alpha, rng=rng) for i in range(n_i)]


# In[6]:


my_model = abm.Model(employer_list, employee_list)
my_alpha, beta, ww, xx, choice, true_opp, obs_opp =     my_model.matching_process(employer_list=employer_list, employee_list=employee_list)
wa = np.array([ee.wa for ee in employee_list])

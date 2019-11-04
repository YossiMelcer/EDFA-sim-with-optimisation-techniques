
# coding: utf-8

# In[47]:


import numpy as np
import matplotlib.pyplot as plt
import scipy.constants
from scipy.integrate import odeint
from functools import partial


# In[81]:


def rate_equations(N,t, sigma_a=2.7e-21,sigma_e=0, phi_p=1,phi_s=1,trans_prob_21=1,trans_prob_32=1):
    N_1 = N[0]
    N_2 =N[1]
    N_3 = N[2]
    dN_1dt = trans_prob_21 * N_2 -(N_1-N_3)*phi_p * sigma_a + (N_2-N_1)*phi_s*sigma_e
    dN_2dt = -trans_prob_21 * N_2 + (-trans_prob_32)*N_3*(N_2-N_1) * phi_s*sigma_e
    dN_3dt = -trans_prob_32 * N_3 + (N_1-N_3)* phi_p*sigma_a
    dNdt = [dN_1dt,dN_2dt,dN_3dt]
    return dNdt


# In[108]:


n_0=[100,1,1]
t=np.linspace(0,10,100)


# In[109]:


N=odeint(rate_equations,n_0,t)


# In[110]:


plt.plot(t,N)


# In[7]:


def n_sp (N_1, N_2, sigma_a, sigma_e):
    return (sigma_e*N_2)/((sigma_e *N_2)-(sigma_a *N_1))


# In[8]:


def P_ASE (n_sp, G, pump_freq, B ):
    return n_sp*(G-1)*scipy.constants.h*pump_freq*B


# In[9]:


def SNR_out (G, I_sig, N_sig_ASE, N_ASE_ASE , N_shot):
    return np.square((G*I_sig)) / (N_sig_ASE + N_ASE_ASE + N_shot)


# In[13]:


def SNR_in(I_sig,B):
    return I_sig/(2*scipy.constants.e*B)


# In[10]:


def Power (signal):
    return np.square(np.abs(signal,dtype=np.float64))


# In[14]:


def F_n(n_sp, G):
    return 2*n_sp*(G-1)/G


# In[15]:


def Gain(sig_out, sig_in):
    return 10*np.log(Power(sig_out)/Power(sig_in))


# In[16]:


Pumping_scheme={980:{'sigma_a': 2.7e-21, 'sigma_e':0}, 1480:{'sigma_a': 1.5e-21, 'sigma_e':0.5e-21}}


# In[87]:


def solve_rate(pump_scheme=980, N_0=[100,100,100], phi_p=1,phi_s=1,
               trans_prob_21=1,trans_prob_32=1, t=np.linspace(0,10,100)):
    
    eq_to_solve = lambda n,t: rate_equations(n, t=t, phi_p=1,phi_s=1,trans_prob_21=1,
                                           trans_prob_32=1, **Pumping_scheme[pump_scheme])
    N=odeint(func=eq_to_solve,y0=N_0, t=t)
    return N


# In[88]:


N= solve_rate()


# In[89]:


plt.plot(N)



# %%


# %%


# %%

# -*- coding: utf-8 -*-
"""
Created on Sun Sep 14 21:53:22 2025

@author: Administrator
"""

import numpy as np
from scipy.integrate import quad
import matplotlib.pyplot as plt

v_f=1/60 # feeding speed
v_d=100/60 # drawing speed
L=0.2 # furnace length
alpha=14300 # temperature profile parameter

# T_max = np.linspace(150, 200, 100) 
# lambda_0 = np.linspace(0.01, 0.1, 100)

def _T(z, T_max): # temperature profile
    return T_max-alpha*(z-L/2)*(z-L/2)

def _eta(z, T_max): # viscosity profile ()
    return np.exp(22493 / ( _T(z, T_max) + 273.15) - 35.287)

def _v(z, T_max): # velocity profile
    _v_denom, _ = quad(lambda g: 1.0 / _eta(g, T_max), 0, L)
    _v_num, _ = quad(lambda g: 1.0 / _eta(g, T_max), 0, z)
    return np.exp(np.log(v_f) + (_v_num/_v_denom) * np.log(v_d / v_f))

def _gamma(z, T_max): #surface tension profile
    return 49.2 - 0.06 * (_T(z, T_max) - 20)

def _lambda(z, T_max, lambda_0): #periodicity profile
    return lambda_0 * 1e-6 * np.sqrt(v_f / _v(z, T_max))

def _tau(z, T_max, lambda_0): #characteristic time
    return _eta(z, T_max) * _lambda(z, T_max, lambda_0) / (3.14 * _gamma(z, T_max))

def _fsh(lambda_0, T_max): #shape factor
    _fsh_func, _ = quad(lambda g: -1.0 / (_tau(g, T_max, lambda_0) * _v(g, T_max)), 0, L)
    return np.exp(_fsh_func)


T_max_list = np.linspace(150, 200, 100)      
lambda_0_list = np.linspace(500, 100, 100)  
fsh_values = np.zeros((len(lambda_0_list), len(T_max_list)))

for i, lam in enumerate(lambda_0_list):
    for j, Tmax in enumerate(T_max_list):
            fsh_values[i, j] = _fsh(lam, Tmax)
            
plt.figure(figsize=(10, 8))
heatmap = plt.imshow(fsh_values.T, aspect='auto', 
                     extent=[lambda_0_list[0], lambda_0_list[-1], 
                             T_max_list[0], T_max_list[-1]],
                     origin='lower', cmap='jet')
plt.colorbar(heatmap, label='Shape Factor (fsh)', ticks=np.linspace(0,1,11))
plt.ylabel('T_max (°C)')
plt.xlabel('λ₀ (μm)')
plt.title('Heatmap of Shape Factor (fsh) vs T_max and λ₀')
plt.tight_layout()
plt.show()

# test
# val=_fsh(100, 150)
# print(val)
# val=_fsh(100, 200)
# print(val)
# val=_fsh(500, 150)
# print(val)
# val=_fsh(500, 200)
# print(val)

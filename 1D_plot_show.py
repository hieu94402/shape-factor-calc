# -*- coding: utf-8 -*-
"""
Created on Wed Sep 10 10:41:33 2025

@author: hieu9
"""

import numpy as np
from scipy.integrate import quad
import matplotlib.pyplot as plt

# define constant
v_f=1/60 # feeding speed
v_d=100/60 # drawing speed
L=0.2 # furnace length
T_max=170 # max temperature, degree
lambda_0=200 #initial periodicity, μm
alpha=14300 # temperature profile parameter

#define variable
z=np.linspace(0, L, 1000)


def _T(z): # temperature profile
    return T_max - alpha * (z-L/2) * (z-L/2)
# _T_vec = np.vectorize(_T)
# _T_values = _T_vec(z)
# plt.figure(figsize=(8,5))
# plt.plot(z, _T_values, marker='o', markersize=1, linestyle='-', color='r')
# plt.title("Temperature vs. position along the fiber")
# plt.xlabel("z coordinate (m)")
# plt.ylabel("T(z) (C degree)")
# plt.grid(True)
# plt.show()

def _eta(z): # viscosity profile ()
    return np.exp(22493 / (_T(z) + 273.15) - 35.287)
def _log_eta(z):
    return np.log(_eta(z))
# _log_eta_vec = np.vectorize(_log_eta)
# _log_eta_values = _log_eta_vec(z)
# plt.figure(figsize=(8,5))
# plt.plot(z, _log_eta_values, marker='o', markersize=1, linestyle='-', color='r')
# plt.title("Viscosity profile (in terms of ln) along the fiber")
# plt.xlabel("z coordinate (m)")
# plt.ylabel("ln η(z)")
# plt.grid(True)
# plt.show()

def _v(z): # velocity profile
    _v_denom, _ = quad(lambda g: 1.0 / _eta(g), 0, L)
    _v_num, _ = quad(lambda g: 1.0 / _eta(g), 0, z)
    return np.exp(np.log(v_f) + (_v_num/_v_denom) * np.log(v_d/v_f))
#_v_vec = np.vectorize(_v)
# _v_values = _v_vec(z)
# plt.figure(figsize=(8,5))
# plt.plot(z, _v_values, marker='o', markersize=1, linestyle='-', color='r')
# plt.title("Velocity profile along the fiber")
# plt.xlabel("z coordinate (m)")
# plt.ylabel("v(z) (mm/s)")
# plt.grid(True)
# plt.show()

def _gamma(z): #surface tension profile
    return 49.2 - 0.06 * (_T(z) - 20)
# _gamma_vec = np.vectorize(_gamma)
# _gamma_values = _gamma_vec(z)
# plt.figure(figsize=(8,5))
# plt.plot(z, _gamma_values, marker='o', markersize=1, linestyle='-', color='r')
# plt.title("Surface tension vs. position along the fiber")
# plt.xlabel("z coordinate (m)")
# plt.ylabel("γ(z) (mN/m)")
# plt.grid(True)
# plt.show()

def _lambda(z): #periodicity profile
    return lambda_0 * np.sqrt(v_f / _v(z))
# _lambda_vec = np.vectorize(_lambda)
# _lambda_values = _lambda_vec(z)
# plt.figure(figsize=(8,5))
# plt.plot(z, _lambda_values, marker='o', markersize=1, linestyle='-', color='r')
# plt.xlim(0, L)
# plt.title("Periodicity along the fiber")
# plt.xlabel("z coordinate (m)")
# plt.ylabel("λ(z) (μm)")
# plt.grid(True)
# plt.show()

def _tau(z): #characteristic time
    return np.log( _eta(z) * 1e-6 * _lambda(z) / (3.14*_gamma(z)) )
# _tau_vec = np.vectorize(_tau)
# _tau_values = _tau_vec(z)
# plt.figure(figsize=(8,5))
# plt.plot(z, _tau_values, marker='o', markersize=1, linestyle='-', color='r')
# plt.title("Characteristic time (in terms of ln) along the fiber")
# plt.xlabel("z coordinate (m)")
# plt.ylabel("ln τ(z)")
# plt.grid(True)
# plt.show()


fsh, _ = quad(lambda g: 1 / (_tau(g) * _v(g)), 0, L)
print(fsh)


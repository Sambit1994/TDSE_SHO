# -*- coding: utf-8 -*-
"""
Created on Sun May 23 09:43:28 2021

@author: HP
"""

##################################################
# Inputs for the module                          #
#--------                                        #
# span_i, span_f : box edges (length of box)     #
# h : step size                                  #
##################################################
##The module creates the 1-D SHO Hamiltonian from the user's input


import numpy as np
import H_module as Hmod
import matplotlib.pyplot as plt


span_i = -6.0
span_f = 6.0

h = 0.1 #0.05


eigvals, eigvecs, X = Hmod.H_mod(span_i, span_f, h)

print ('1st 5 eigen values are:', '%7.5f \t %7.5f \t %7.5f \t %7.5f \t %7.5f' % (eigvals[0], eigvals[1], eigvals[2], eigvals[3], eigvals[4]))
   
figg, ax = plt.subplots(figsize=(6, 6))
plt.xlabel("x",fontsize=15)
plt.ylabel("\u03C8 (x)",fontsize=15)
#ax.plot(H_set1, rho_set1)
ax.plot(X, eigvecs[:,0], label="state-1") 
ax.plot(X, eigvecs[:,1], label="state-2")
ax.plot(X, eigvecs[:,2], label="state-3")
ax.plot(X, eigvecs[:,3], label="state-4")
ax.legend(loc='upper right')
plt.show()
figg.savefig('Eigen_Vectors.svg', dpi=600)

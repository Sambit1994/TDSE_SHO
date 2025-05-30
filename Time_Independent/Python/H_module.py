# -*- coding: utf-8 -*-
"""
Created on Sun May 23 09:43:23 2021

@author: HP
"""

######################################################################
# Inputs                                                             #
#--------                                                            #
# span_i, span_f : box edges (length of box)                         #
# h : step size                                                      #
# KE -- Five-point central difference formula for second derivatives #
# VE (w/o bump) -- x**2/2                                            #
######################################################################

##The module creates the 1-D SHO Hamiltonian from the user's input

import numpy as np
import math
import matplotlib.pyplot as plt
from scipy.linalg import eigh

def H_mod(span_i,span_f,h):

   x_val = span_i - h
   
   N_it = ((span_f-span_i)/h) + 1
   
   N_end   = int(N_it)
   
   V_const = 1.0/(24.0 * (h**2))
   
   V1 = 1.0 * V_const
   V2 = -16.0 * V_const
   V3 = 30.0 * V_const 
   V4 = -16.0 * V_const
   V5 = 1.0 * V_const
   
   
   S = np.zeros((N_end, N_end)) 
   #print (S)
   
   N_start = 0
   
   X = []
   for i in range(0, (N_end)):
       x_val = x_val + h
       X = np.append(X, x_val)
       #print (x_val)
   
       f0 = 2*(np.exp(-10*(x_val**2)))  
       f1 = 3*(np.exp(-15*(x_val**2)))
       f2 = 2*(np.exp(-5*(x_val**2)))
       f3 = 2*(np.exp(-30*(x_val**2)))
       f4 = 10*(np.exp(-30*(x_val**2)))
       f5 = 10*(np.exp(-5*(x_val**2)))
       f5_2 = 10*(np.exp(-5*((x_val-1)**2))) + 10*(np.exp(-5*((x_val+1)**2)))
       f5_3 = 10*(np.exp(-5*((x_val-1)**2))) + 10*(np.exp(-5*((x_val)**2))) + 10*(np.exp(-5*((x_val+1)**2)))
       f5_4 = 10*(np.exp(-5*((x_val-1)**2)))
       f6 = 30*(np.exp(-5*(x_val**2)))

       V_pot = ((x_val**2)/2) + f5
       
       j = i-2
       k = i-1
       l = i+1
       m = i+2
       
       
       if ((k == (N_start-1)) and (j == (N_start-2))):
           S[i,i] = V3 + V_pot
           S[i,l] = V4
           S[i,m] = V5        
       elif ((k == (N_start)) and (j == (N_start-1))):
           S[i,k] = V2
           S[i,i] = V3 + V_pot
           S[i,l] = V4
           S[i,m] = V5
       elif ((l == (N_end)) and (m == (N_end+1))):
           if (m < N_end):
               S[i,j] = V1
               S[i,k] = V2
               S[i,i] = V3 + V_pot
               S[i,l] = V4
           else:
               S[i,j] = V1
               S[i,k] = V2
               S[i,i] = V3 + V_pot            
       elif ((l == (N_end+1)) and (m == (N_end+2))):
           S[i,j] = V1
           S[i,k] = V2
           S[i,i] = V3 + V_pot
       else:
           if (m < N_end):
               S[i,j] = V1
               S[i,k] = V2
               S[i,i] = V3 + V_pot
               S[i,l] = V4
               S[i,m] = V5
           else:
               S[i,j] = V1
               S[i,k] = V2
               S[i,i] = V3 + V_pot
               S[i,l] = V4            
   
   eigvals, eigvecs = eigh(S)
   
   return (eigvals, eigvecs, X)


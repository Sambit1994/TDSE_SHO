# -*- coding: utf-8 -*-
"""
Created on Sat May 29 08:11:53 2021

@author: HP
"""

######################################################################
# Variables                                                          #
#----------                                                          #
# span_i, span_f : box edges (length of box)                         #
# h : step size                                                      #
# Nsegments : No. of time steps                                      #
# Omega : Larger Omega --> Smaller time period                       #
# SHO Hamiltonian uses                                               #
#----------------------                                              #
# KE -- Five-point central difference formula for second derivatives #
# VE (w/o bump) -- x**2/2                                            #
######################################################################


import numpy as np
import math
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import matplotlib.animation as animation
from scipy.linalg import eigh
#from IPython import display
#import subprocess

span_i = -6.0
span_f = 6.0

h = 0.1 #0.05
x_val = span_i - h

N_it = ((span_f-span_i)/h) + 1

N_start = 1
N_end   = int(N_it)

KE_const = 1.0/(24.0 * (h**2))

KE1 = 1.0 * KE_const
KE2 = -16.0 * KE_const 
KE3 = 30.0 * KE_const 
KE4 = -16.0 * KE_const
KE5 = 1.0 * KE_const

#### Function Creating Hamiltonian matrix S and Potential ####

S = np.zeros((N_end, N_end)) 

N_start = 0

Omega = 2
t_last = np.pi/Omega
#print ('t_last',t_last)

Nsegments = 51
h_t = t_last/(Nsegments -1)

t = np.arange(0, (t_last+h_t), h_t)
#print ('t',len(t), t)

def solve_H(time, x_val):
    X = []
    V_POTEN = []
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
        f_t = (np.sin(Omega*time))*f5_3
    
        V_pot = ((x_val**2)/2) + f_t 
        
        #if (time > t_last):
        #    #print(time)
        #    V_pot = ((x_val**2)/2)
        #else:
        #    #print ('within range', time)            
        #    V_pot = ((x_val**2)/2) + f_t
        
        V_POTEN = np.append(V_POTEN, V_pot)
        j = i-2
        k = i-1
        l = i+1
        m = i+2
          
        if ((k == (N_start-1)) and (j == (N_start-2))):
            S[i,i] = KE3 + V_pot
            S[i,l] = KE4
            S[i,m] = KE5        
        elif ((k == (N_start)) and (j == (N_start-1))):
            S[i,k] = KE2
            S[i,i] = KE3 + V_pot
            S[i,l] = KE4
            S[i,m] = KE5
        elif ((l == (N_end)) and (m == (N_end+1))):
            if (m < N_end):
                S[i,j] = KE1
                S[i,k] = KE2
                S[i,i] = KE3 + V_pot
                S[i,l] = KE4
            else:
                S[i,j] = KE1
                S[i,k] = KE2
                S[i,i] = KE3 + V_pot            
        elif ((l == (N_end+1)) and (m == (N_end+2))):
            S[i,j] = KE1
            S[i,k] = KE2
            S[i,i] = KE3 + V_pot
        else:
            if (m < N_end):
                S[i,j] = KE1
                S[i,k] = KE2
                S[i,i] = KE3 + V_pot
                S[i,l] = KE4
                S[i,m] = KE5
            else:
                S[i,j] = KE1
                S[i,k] = KE2
                S[i,i] = KE3 + V_pot
                S[i,l] = KE4            

    return (X, S, V_POTEN)        
        
#print (S)      
#### Function Creating Hamiltonian matrix S and Potential #### --End


#### Check Point1 --- Validating the calcuated eigen values (at t=0; w/o bump) are correct ####

X, S, V_POTEN0 = solve_H(t[0], x_val)
eigvals, eigvecs = eigh(S)

eigvecs1 = eigvecs.transpose()
pdt = np.matmul(eigvecs, eigvecs1)

CHK1 = np.matmul((eigvecs.transpose()), S)
CHK = np.matmul(CHK1, eigvecs) ## Checking the eigen values in matrix format <Psi|H|Psi>
#print (CHK)
#print ('CHK')
#print (CHK[0,0], CHK[10,10], CHK[12,12], CHK[21,21])

print ('1st 5 eigen values (at t=0) are:', '%7.5f \t %7.5f \t %7.5f \t %7.5f \t %7.5f' % (eigvals[0], eigvals[1], eigvals[2], eigvals[3], eigvals[4]))

#### Check Point1 #### ---End


#### Function for generating time-dependent Hamiltonian matrix and Potential ####

#print ('t',t)
n_state = 0

eigvecs_oldd = eigvecs.copy()
eigvecs_old = eigvecs_oldd[:,n_state]

def solve_H_T(t_new, eigvecs_old):    
    X, S, V_POTEN_new = solve_H(t_new, x_val)   
    eigvals_new, eigvecs_new = eigh(S)
    eigvecs_new_T = eigvecs_new.transpose()
    #print('1st 5 eigen values at t=%7.5f are: %7.5f, %7.5f, %7.5f, %7.5f, %7.5f' % (t_new, eigvals_new[0], eigvals_new[1], eigvals_new[2], eigvals_new[3], eigvals_new[4]))
    Psi_total = 0
    for k in range(len(eigvecs_new[0])):
        
        pd1 = np.matmul(eigvecs_new[:,k].transpose(), eigvecs_old)            
        
        D = np.dot(pd1, eigvecs_new[:,k])
        
        value = eigvals_new[k]*h_t
        
        vl = complex(np.cos(value), -(np.sin(value)))
        
        Psi_t = np.dot(vl, D)
        #print ('Psi-t',Psi_t)
        Psi_total = Psi_total + Psi_t
        #print (Psi_total)

    eigvecs_old = Psi_total.copy()        
    return Psi_total, V_POTEN_new

#### Function for generating time-dependent Hamiltonian matrix and Potential #### --End



#### Check Point2 -- plotting time dependent potential (Fig-1:ax) ####

figg, ax = plt.subplots(figsize=(6, 6))
plt.xlabel("x",fontsize=15)
plt.ylabel("V(x,t)",fontsize=15)

ax.plot(X, V_POTEN0, label='t0') 

eigvecs_old_set = np.zeros([len(t), len(eigvecs_oldd[:,0])],dtype = 'complex_')
eigvecs_old_set[0] = eigvecs_oldd[:,n_state]

for l in range(len(t)-1):
    m = l+1
    m_i = "{:.3f}".format(m)
    #print (m)
    Psi_total, V_POTENn = solve_H_T(t[m], eigvecs_old)
    eigvecs_old = Psi_total.copy()
    eigvecs_old_set[m] = eigvecs_old.copy()
    ax.plot(X, V_POTENn, label=f"t{m}")

ax.set_title("Time evolution of Potential")
#ax.legend(loc='upper right')

#### Check Point2 #### --End



#### plotting time dependent potential Animation (Fig-2:bx) ####

x_data = X.copy()
y_data = V_POTEN0

fig, bx = plt.subplots(figsize=(6, 5))
plt.xlabel("x",fontsize=15)
plt.ylabel("V(x,t)",fontsize=15)

bx.set_title("Time evolution of Potential (Animation)") 

line, = bx.plot(x_data, y_data, linewidth=4)

n_frames = np.delete(t, np.s_[0], axis=0)

eigvecs_old = eigvecs_oldd[:,n_state]

def animation_frame(l):
    index, = np.where(n_frames==l)
    ind = int(index)
    Psi_total, V_POTENn = solve_H_T(l, eigvecs_old_set[ind])
    x_data = X.copy()
    y_data = V_POTENn
    line.set_xdata(x_data)
    #line.set_ydata(abs(y_data))
    line.set_ydata(y_data)
    return line, 

anim = FuncAnimation(fig, func=animation_frame, frames=n_frames, interval=100, blit=True)

# use 'matplotlib qt' command in the console window to run the animation in spyder

#### plotting time dependent potential Animation (Fig-2:bx) #### --End

writergif = animation.PillowWriter(fps=100) 
anim.save('Animation_pot_function.gif', writer=writergif)

plt.show()
figg.savefig('Potential_function.svg', dpi=600)




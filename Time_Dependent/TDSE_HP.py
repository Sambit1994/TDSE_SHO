# -*- coding: utf-8 -*-
"""
Created on Mon May 31 16:54:36 2021

@author: HP
"""

##################################################
# Variables                                      #
#----------                                      #
# span_i, span_f : box edges (length of box)     # 
# h : step size                                  #
# Nsegments : No. of time steps                  #
# Omega : Larger Omega --> Smaller time period   #
##################################################


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

h = 0.1 #0.05 #Step size
x_val = span_i - h

#No of iterations
N_it = ((span_f-span_i)/h) + 1

N_start = 1
N_end   = int(N_it)


KE_const = 1.0/(24.0 * (h**2))

KE1 = 1.0 * KE_const
KE2 = -16.0 * KE_const #-2.0 #-16.0
KE3 = 30.0 * KE_const #1.0 #30.0
KE4 = -16.0 * KE_const
KE5 = 1.0 * KE_const


#### Function Creating Hamiltonian matrix S ####

S = np.zeros((N_end, N_end)) 
#print (S)

N_start = 0

Omega = 2#20#0.2
t_last = np.pi/Omega
#t_last = (2*np.pi)/Omega
print (t_last)

Nsegments = 51
h_t = t_last/(Nsegments -1)
t = np.arange(0, (t_last+h_t), h_t)
#print ('t',len(t), t)

def solve_H(time, x_val):
    X = []
    for i in range(0, (N_end)):
        x_val = x_val + h
        X = np.append(X, x_val)
        #print (x_val)
    
        #f0 = 1.0
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
        f_t = (np.sin(Omega*time))*f5
    
        #V_pot = ((x_val**2)/2) + f_t 

        if (time > t_last):
            #print(time)
            V_pot = ((x_val**2)/2)
        else:
            #print ('within range', time)            
            V_pot = ((x_val**2)/2) + f_t

        
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

    return (X, S)        
        
#### Function Creating Hamiltonian matrix S #### --End


#### Check Point1 --- Validating the calcuated eigen values (at t=0; w/o bump) are correct ####
X, S = solve_H(t[0], x_val)
eigvals, eigvecs = eigh(S)

eigvecs1 = eigvecs.transpose()
pdt = np.matmul(eigvecs, eigvecs1)

#CHK1 = np.matmul((eigvecs.transpose()), S)
#CHK = np.matmul(CHK1, eigvecs) ## Checking the eigen values in matrix format <Psi|H|Psi>
#print (CHK)
#print ('CHK')
#print (CHK[0,0], CHK[1,1], CHK[10,10], CHK[12,12], CHK[21,21])

print ('1st 5 eigen values (at t=0) are:', '%7.5f \t %7.5f \t %7.5f \t %7.5f \t %7.5f' % (eigvals[0], eigvals[1], eigvals[2], eigvals[3], eigvals[4]))

#### Check Point1 #### ---End



#### Function for generating time-dependent Hamiltonian matrix ####

#print ('t',t)
n_state = 0

eigvecs_oldd = eigvecs.copy()
eigvecs_old = eigvecs_oldd[:,n_state]

def solve_H_T(t_new, eigvecs_old):    
    X, S = solve_H(t_new, x_val)   
    eigvals_new, eigvecs_new = eigh(S)
    eigvecs_new_T = eigvecs_new.transpose()
    print('1st 5 eigen values at t=%7.5f are: %7.5f, %7.5f, %7.5f, %7.5f, %7.5f' %
      (t_new, eigvals_new[0], eigvals_new[1], eigvals_new[2], eigvals_new[3], eigvals_new[4]))
    
    Psi_total = 0
    pd1_set = []
    for k in range(len(eigvecs_new[0])):
        
        pd1 = np.matmul(np.conj(eigvecs_new[:,k]).transpose(), eigvecs_old)            
        pd1_set = np.append(pd1_set, pd1)

        D = np.dot(pd1, eigvecs_new[:,k])
        
        value = eigvals_new[k]*h_t
        
        vl = complex(np.cos(value), -(np.sin(value)))
        
        Psi_t = np.dot(vl, D)

        #print ('Psi-t',Psi_t)
        Psi_total = Psi_total + Psi_t  ##Time evolution of wave function

    eigvecs_old = Psi_total.copy() 
    return Psi_total, eigvecs_new, pd1_set

#### Function for time-dependent Hamiltonian matrix #### --End


#### Check Point2 -- plotting eigen function at t=0 (Fig-1:ax) ####

figg, ax = plt.subplots(figsize=(6, 6))
plt.xlabel("x",fontsize=15)
plt.ylabel("\u03C8 (x)",fontsize=15)

ax.plot(X, np.real(eigvecs[:,n_state]), label="state-1")  ##Use this
ax.plot(X, np.real(eigvecs[:,(n_state+1)]), label="state-2") 
ax.plot(X, np.real(eigvecs[:,(n_state+2)]), label="state-3") 
ax.plot(X, np.real(eigvecs[:,(n_state+3)]), label="state-4") 

ax.set_title("Eigenvectors (States 1-4) at t=0")
ax.legend(loc='upper right')

#### Check Point2 #### --End


#### Determining time evolution of wave function and plotting ####

#Plotting the wavepacket or probability densities (commented out) over time (Fig-2:ax2)

figg2, ax2 = plt.subplots(figsize=(6, 6))
plt.xlabel("x",fontsize=15)
plt.ylabel("|\u03C8 (x)|",fontsize=15) ##to plot |Psi>
#plt.ylabel("|\u03C8 (x)|2",fontsize=15) ##to plot |Psi>**2

ax2.plot(X, np.real(eigvecs[:,n_state]))  ##to plot |Psi>
#ax2.plot(X, np.real(eigvecs[:,n_state]*np.conj(eigvecs[:,n_state])))  ##to plot |Psi>**2 

eigvecs_old_set = np.zeros([len(t), len(eigvecs_oldd[:,0])],dtype = 'complex_')
eigvecs_old_set[0] = eigvecs_oldd[:,n_state]

#Plotting coefficients at each time step (Fig-3:dx)

## coefficient at each time step <i(t+1)|psi(t)>
state_no = np.arange(1, (len(eigvecs_oldd[:,0])+1), 1)
#print (state_no)
fig_cont_each, dx = plt.subplots(figsize=(6, 6))
plt.xlabel("state no",fontsize=15)
plt.ylabel("|C|$^2$",fontsize=15)
dx = plt.axes(xlim=(0, 20))

for l in range(len(t)-1):
    m = l+1
    Psi_total, eigvecs_new, pd1_set = solve_H_T(t[m], eigvecs_old)
    eigvecs_old = Psi_total.copy()
    eigvecs_old_set[m] = eigvecs_old.copy()
    
    ax2.plot(X, np.real(Psi_total)) ##to plot |Psi>
    #ax2.plot(X, np.real(Psi_total*np.conj(Psi_total))) ##to plot |Psi>**2
    
    #print (Psi_total*np.conj(Psi_total)) ##coefficients
    dx.plot(state_no, np.real(pd1_set*np.conj(pd1_set)), marker='o', linestyle='dashed')

ax2.set_title("Time evolution of wavepacket") ##to plot |Psi>
#ax2.set_title("Time evolution of Probability density") ##to plot |Psi>**2
#ax2.legend(loc='upper right')

dx.set_title("Contribution of eigenstates at each time step")
#dx.legend(loc='upper right')


## soving Hamiltonian over the time and plotting time evolution of eigen funtion:
fig_each_H, ex = plt.subplots(figsize=(6, 6))
plt.xlabel("x",fontsize=15)
plt.ylabel("|\u03C8 (x)|2",fontsize=15)
for o in range(len(t)):
    #print (t[o])
    X, S = solve_H(t[o], x_val)
    eigvals_H_new, eigvecs_H_new = eigh(S)
    #ex.plot(X, np.real(eigvecs_H_new[:,0]), label=f"t = {t[o]:.5f}")   #State-1 
    #ex.plot(X, np.real(eigvecs_H_new[:,1]), label=f"t = {t[o]:.5f}")   #State-2
    ex.plot(X, np.real(eigvecs_H_new[:,2]), label=f"t = {t[o]:.5f}")    #State-3

ex.set_title("Time evolution of eigenvector (State-3)")
#ex.legend(loc='upper right')

#print ('New eigen values', eigvals_new)

#### Determining time evolution of wave function and plotting #### --End 



#### Expectation value of the final wavefunction and internal checking ####

eigvecs_old = eigvecs_old_set[-1].copy()
X1, S1 = solve_H(t[0], x_val)   
eigvals_last, eigvecs_last = eigh(S1)
M1 = np.matmul(eigvecs_last[:,0].transpose(), S1)
M2 = np.matmul(M1, eigvecs_last[:,0])
#print ('M2',M2)
#print ('chk')
pd1_set1 = []
for k in range(len(eigvecs_new[0])):
    pd11 = np.matmul(np.conj(eigvecs_last[:,k]).transpose(), eigvecs_old)            
    pd1_set1 = np.append(pd1_set1, pd11)

#print ('pd1_set1')
#print (pd1_set1)
#print ()
#print (pd1_set1.real)
#print ('pd1_set1',(np.real(pd1_set1*np.conj(pd1_set1))), sum(np.real(pd1_set1*np.conj(pd1_set1))))
#print ('pd1_set1',(np.matmul(pd1_set1,np.conj(pd1_set1))))

AA1 = np.real(pd1_set1*np.conj(pd1_set1))
#print ()

eigvecs_old = eigvecs_old_set[-1].copy()
X2, S2 = solve_H(t[-1], x_val) 
#X2, S2 = solve_H(t[0], x_val) ## Both will work as both corresponds to the normal Hamiltonian
evals, evecs = eigh(S2)
#print (evals)
AA2 = evals
#print ()
M11 = np.matmul(eigvecs_old.transpose(), S2)
M22 = np.matmul(M11, np.conj(eigvecs_old))
#print ('M22',M22)
#print (M22*np.conj(M22))
#print ('check')
#print (np.matmul(AA1,AA2))
#print ()

    
#print (len(eigvecs_old_set))
X1, S1 = solve_H(t[0], x_val)   
eigvals_HO, eigvecs_HO = eigh(S1)
for i in range(len(eigvecs_old_set)):
    eigvecs_old = eigvecs_old_set[i]
    pd1_set1 = []
    for k in range(len(eigvecs_HO[0])):
        #pd11 = np.matmul(eigvecs_HO[:,k].transpose(), eigvecs_old)
        pd11 = np.matmul(np.conj(eigvecs_HO[:,k]).transpose(), eigvecs_old)            
        pd1_set1 = np.append(pd1_set1, pd11)
    #print ('pd1_set1')
    #print (pd1_set1)
    coeffcients_sq = np.real(pd1_set1*np.conj(pd1_set1))
    vals = eigvals_HO.copy()
    result = np.matmul(coeffcients_sq,vals)
    #print ('coeff', coeffcients_sq[2])
    #print (result)

print ()    
print ('Final expectation value of the resultant wave function is:', '%7.5f' % (result))

#### Expectation value of the final wavefunction and internal checking #### --End 




#### Plotting contribution of eigenstates to the final wave function ####

state_no = np.arange(1, (len(coeffcients_sq)+1), 1)
#print (state_no)
fig_satte, cx = plt.subplots(figsize=(6, 6))
plt.xlabel("state no",fontsize=15)
#plt.ylabel("\u03C8 (x)",fontsize=15)
plt.ylabel("|C|$^2$",fontsize=15)
cx = plt.axes(xlim=(0, 20))
cx.plot(state_no, coeffcients_sq, color='green', marker='o', linestyle='dashed')

cx.set_title("Contribution of eigenstates to the final wave function")

#### Plotting contribution of eigenstates to the final wave function #### --End



#### Plotting time evolution of wavepacket or probability density (commented out) as animation ####

x_data = X.copy()
y_data = eigvecs[:,n_state].copy()
y_data = np.real(y_data) ##To plot |Psi>
#y_data = np.real(y_data*(np.conj(y_data))) ##To plot |Psi>**2

fig, bx = plt.subplots(figsize=(6, 5))
plt.xlabel("x",fontsize=15)

plt.ylabel("|\u03C8 (x)|",fontsize=15) ##To plot |Psi>
#plt.ylabel("|\u03C8 (x)|2",fontsize=15) ##To plot |Psi>**2

bx = plt.axes(ylim=(-0.2, 0.25)) ##To plot |Psi>
bx.set_title("Time evolution of wave packet (Animation)") ##To plot |Psi>
#bx.set_title("Time evolution of probability density (Animation)") ##To plot |Psi>**2
#bx = plt.axes(ylim=(-0.02, 0.1)) ##To plot |Psi>**2

line, = bx.plot(x_data, y_data, linewidth=4) 

n_frames = np.delete(t, np.s_[0], axis=0)

eigvecs_old = eigvecs_oldd[:,n_state]

def animation_frame(l):
    index, = np.where(n_frames==l)
    ind = int(index)
    x_data = X.copy()
    y_data = eigvecs_old_set[ind+1]
    y_data = np.real(y_data) ##To plot |Psi>
    #y_data = np.real(y_data*(np.conj(y_data))) ##To plot |Psi>**2 
    line.set_xdata(x_data)
    #line.set_ydata(abs(y_data))
    line.set_ydata(y_data)
    return line, 

#anim = FuncAnimation(fig, func=animation_frame, frames=np.arange(0, 10, 0.1), interval=10)
anim = FuncAnimation(fig, func=animation_frame, frames=n_frames, interval=100, blit=True)
#anim = FuncAnimation(fig, func=animation_frame, frames=n_frames, interval=100, repeat=False)
#anim = FuncAnimation(fig, func=animation_frame, interval=10)

# use 'matplotlib qt' command in the console window to run the animation in spyder

#### Plotting time evolution of wavepacket or probability density (commented out) as animation #### --End


writergif = animation.PillowWriter(fps=100) 
anim.save('animation.gif', writer=writergif)

plt.show()
figg.savefig('Eigenvector_t0.svg', dpi=600)
figg2.savefig('Time_evolution_wavepacket.svg', dpi=600)
fig_cont_each.savefig('Eigenstate_cont_vs_time.svg', dpi=600)
fig_each_H.savefig('Eigenvector_time_evolv.svg', dpi=600)
fig_satte.savefig('Eigenstate_cont_final.svg', dpi=600)


The TDSE_SHO.py routine solves the time dependent Schrödinger equation and determines the evolution of eigenvalues and eigenvectors for a 1-D SHO. A bump function ($10e^{-5x^{2}}$) is added to the potential.

---

The **Time Dependent Schrödinger Equation** takes the form

$$\Psi(x,t) = \sum_{n=1}^N C_n(0) e^{-iE_nt/\hbar} \psi_n(x)$$, where contribution from various states can be determined as:

$$C_n = \braket{\psi_n(x)}{\Psi(x,t)}$$

The TDSE_SHO.py adds this on top of the routine presented in the `Time-Independent` formalism.

The time-dependent Hamiltonian can be presented as:

$$H(t) = H_\omega(x) + V(x,t)$$

$H_\omega(x)$ is the time-independent harmonic oscillator hamiltonian and the $V(x,t)$ contains the bump potential with evolving time, i.e.,

$$V(x,t) = \sin(\Omega t)V(x), ~~ 0 \leq t \leq \dfrac{\pi}{\Omega}$$

so, $V(x,t)$ is non-zero only between $t=0~-~\dfrac{\pi}{\Omega}$. Starting in the ground state of the harmonic oscillator at $t_0$ and using a time grid $t_0, t_1, t_2, ...t_n, t_{n+1}, ...$, we will get:

$$\Psi(t_{n+1}) = \sum_i e^{-iE_i(t_{n+1})\Delta t/\hbar} \ket{i_{t+1}}\braket{i_{t+1}}{\Psi(t_n)}$$ --- the wavepacket

---

A series of output files is plotted

**Eigenvector_t0.svg** : Eigenvector at t=0

**Time_evolution_wavepacket.svg** : The evolution of the wavepacket across timesteps

**Eigenstate_cont_vs_time.svg** : Contribution of various eigenstates to the resulting wavepacket across timesteps

**Eigenvector_time_evolv.svg** : Evolution of eigenvector (state-3 as example) across timesteps

**Eigenstate_cont_final.svg** : Contribution of various eigenstates to the resulting wavepacket at the end of simulation

**animation.gif** : A gif image presenting the dynamic evolution of the wavepacket

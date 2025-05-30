Solution of the time dependent Schrödinger equation for a 1D simple harmonic oscillator problem
------------------------------------------------------------------------------------------------
The scripts solve the Schrödinger equation and determine the eigenvalues and eigenvectors for a 1-D simple harmonic oscillator (SHO) problem.

# Quantum Harmonic Oscillator: Hamiltonian and Eigenvalues

This project focuses on computing the **eigenvalues (energies)** and **eigenvectors** of the **Hamiltonian operator** \( H \), which for a one-dimensional system is given by:

$$
H = -\frac{\hbar^2}{2m} \frac{\partial^2}{\partial x^2} + V(x)
$$

- The first term represents the **kinetic energy**.
- The second term represents the **potential energy** \( V(x) \) of the system.

---

## The One-Dimensional Harmonic Oscillator

The **simple harmonic oscillator** serves as a foundational model to describe molecular vibrations. Its importance stems from the fact that **near any local potential minimum**, the system behaves approximately like a harmonic oscillator.

### Taylor Expansion of Potential

Expanding \( V(x) \) about a minimum at \( x = x_0 \):

$$
V(x) = V(x_0) + V'(x - x_0) + \frac{1}{2}V''(x - x_0)^2 + \dots
$$

Ignoring constant and linear terms (assuming \( x_0 = 0 \)):

$$
V(x) \approx \frac{1}{2} m \omega^2 x^2
$$

---

## Final Form of the Hamiltonian

Substituting into the original Hamiltonian:

$$
H = -\frac{\hbar^2}{2m} \frac{\partial^2}{\partial x^2} + \frac{1}{2} m \omega^2 x^2
$$

Introduce a **dimensionless variable**:

$$
z = x \sqrt{\frac{m \omega}{\hbar}}, \quad \text{so that} \quad \frac{\partial^2}{\partial x^2} = \frac{m \omega}{\hbar} \frac{\partial^2}{\partial z^2}
$$

Then the Hamiltonian becomes:

$$
H = \frac{\hbar \omega}{2} \left( -\frac{\partial^2}{\partial z^2} + z^2 \right)
$$

---

## The Eigenvalue Problem

The Schrödinger equation becomes:

$$
\frac{1}{2} \left( -\frac{\partial^2}{\partial z^2} + z^2 \right) \Psi = \frac{E}{\hbar \omega} \Psi
$$

Thus, the **energy eigenvalues** are expressed in **units of \( \hbar \omega \)**.



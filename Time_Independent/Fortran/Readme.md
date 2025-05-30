The H_matrix_finite_diff.f90 routine generates a [n*n] Hamiltonian matrix for a box length ($z_{min}=-10$ to $z_{max}=-10$), where $n = \dfrac{zmax-zmin}{h} + 1 $. $h$ is the step size (0.1). The output Hamiltonian matrix is written to the **mat.dat** file.

The lapak_diagonalisation.f90 routine uses the LAPACK routine $dsyev$ to find the eigenvalues and eigenvectors from **mat.dat** and writes them to the **Diag_mat.dat** file

The first five eigenvectors are collected in the **Eigenvectors.dat** file and plotted (gnuplot plot_eigenvector.plt). The first column is the box length (-10 to 10), spanned with h=0.1.

The Inv_pow_Ite.f90 routine performs the **Inverse Power Iteration with Shift** approach to converge to the nearest eigenvalue. A default shift of 1.3 is used as an example. The output is collected in the **Dia_IPI.dat** file.

--


#Inverse Power Iteration with shift
---

The Power Iteration method provides the highest eigenvalue of a matrix $\bold{A}$. Starting from a guess vector $\bold{Y}_1$,

$\bold{Y}_1 = \sum_n C_n \bold{X}_n$ (since, eigenvectors (still unknown) to $\bold{A}$ form a complete set.)

Applying $\bold{A}$ from left, we get

$\bold{Y}_2 = \bold{A} \bold{Y}_1 = \sum_n C_n \lambda_n \bold{X}_n$

$\bold{Y}_3 = \bold{A} \bold{Y}_2 = \sum_n C_n \lambda_n^2 \bold{X}_n$

In general: $\bold{Y}_{i+1} = \bold{A} \bold{Y}_i = \sum_n C_n \lambda_n^i \bold{X}_n$

After sufficiently many applications of $\bold{A}$ the largest absolute
eigenvalue will dominate completely

$\dfrac{\bold{Y}_i^\dagger \bold{A} \bold{Y}_i}{\bold{Y}_i^\dagger \bold{Y}_i} \Rightarrow \max |\lambda_n|$

For normalised vectors, $\bold{Y}_i^\dagger \bold{A} \bold{Y}_i \Rightarrow \max |\lambda_n|$

---

In contrast to the power iteration, the inverse power iteration calculates the lowest eigenvalue of matrix $\bold{A}$ as:

```math
\bold{A}\bold{X}_n = \lambda_n \bold{X}_n  \Rightarrow \lambda^{-1}\bold{X}_n = \bold{A}^{-1} \bold{X}_n

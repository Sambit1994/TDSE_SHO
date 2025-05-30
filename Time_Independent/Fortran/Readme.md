The H_matrix_finite_diff.f90 routine generates a [n*n] Hamiltonian matrix for a box length ($z_{min}=-10$ to $z_{max}=-10$), where $n = \dfrac{zmax-zmin}{h} + 1 $. $h$ is the step size (0.1). The output Hamiltonian matrix is written to the **mat.dat** file.

The lapak_diagonalisation.f90 routine uses the LAPACK routine $dsyev$ to find the eigenvalues and eigenvectors from **mat.dat** and writes them to the **Diag_mat.dat** file

The first five eigenvectors are collected in the **Eigenvectors.dat** file and plotted (gnuplot plot_eigenvector.plt). The first column is the box length (-10 to 10), spanned with h=0.1.

The Inv_pow_Ite.f90 routine performs the **Inverse Power Iteration with Shift** approach to converge to the nearest eigenvalue. A default shift of 1.3 is used as an example. The output is collected in the **Dia_IPI.dat** file.

--


#Inverse Power Iteration with shift
---

The Power Iteration method provides the highest eigenvalue of a matrix **A**. Starting from a guess vector **Y$_1$**, 

$**Y**_1 = \sum_n C_n **X**_n$ (since, eigenvectors (still unknown) to **A** form a complete set.)

In general: $**Y**_{i+1}$ =  $**A** **Y**_i = \sum_n C_n \lambda_n^i **X**_n$ 

After sufficiently many applications of **A**, the largest absolute eigenvalue will dominate completely 

$\dfrac{**Y**_i^\dagger **A**  **Y**_i}{**Y**_i^\dagger **Y**_i} \Rightarrow max |\lambda_n| $

for normalised vectors, $**Y**_i^\dagger **A** **Y**_i \Rightarrow max |\lambda_n| $

In contrast to the power iteration, the **inverse power iteration** calculates the lowest eigenvalue of matrix **A** as: 

$**A** **X**_n = \lambda_n **X**_n  \Rightarrow \lambda^{-1}**X**_n = **A**^{-1} **X**_n$

where, power iteration of $**A**^{-1}$ gives the lowest eigenvalue. However, instead of determining the $**A**^{-1}$ the following matrix equation is solved 

$**A** **Y**_{i+1} = **Y**_i$\\

Here, one starts from the guess vector **Y**$_i$ and solves the matrix equation to get to the next vector $**Y**_{i+1}$, which further goes through the same procedure. After numerous iterations.

$**Y**_i^\dagger **A**^{-1} **Y**_i \Rightarrow \dfrac{1}{min |\lambda_n|} $; for normalised vectors 

The inverse power iteration with shift is utilised to obtain the smallest eigenvalue of a shifted matrix **A**, such that 

$**A** \Rightarrow **A** - \eta **I**$, where $\eta$ is the scalar constant close to the nearest eigenvalue.

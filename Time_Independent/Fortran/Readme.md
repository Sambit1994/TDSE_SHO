The H_matrix_finite_diff.f90 routine generates a [n*n] Hamiltonian matrix for a box length ($z_{min}=-10$ to $z_{max}=-10$), where $n = \dfrac{zmax-zmin}{h} + 1 $. $h$ is the step size (0.1). The output Hamiltonian matrix is written to the **mat.dat** file. The routine also has the possibility to add a `bump` to the potential and analyze the results. An example is in the `Example_w_bump_potential` section, based on the bump function $10 e^{-5x^{2}}$

The lapak_diagonalisation.f90 routine uses the LAPACK routine $dsyev$ to find the eigenvalues and eigenvectors from **mat.dat** and writes them to the **Diag_mat.dat** file

The first five eigenvectors are collected in the **Eigenvectors.dat** file and plotted (gnuplot plot_eigenvector.plt). The first column is the box length (-10 to 10), spanned with h=0.1.

The Inv_pow_Ite.f90 routine performs the **Inverse Power Iteration with Shift** approach to converge to the nearest eigenvalue. A default shift of 1.3 is used as an example. The output is collected in the **Dia_IPI.dat** file.

---

#Compilation
`gfortran H_matrix_finite_diff.f90 -o H_matrix.x`

`gfortran lapak_diagonalisation.f90 -o lpk_dia.x -llapack -lblas`

`gfortran Inv_pow_Ite.f90 -o Inv_PI.x -llapack -lblas`

---

#Inverse Power Iteration with shift
---

The Power Iteration method provides the highest eigenvalue of a matrix $A$. Starting from a guess vector $Y_1$, 

$Y_1 = \sum_n C_n X_n$ (since, eigenvectors (still unknown) to $A$ form a complete set.) 

Applying $A$ from left, we get 

$Y_2$ =  $A Y_1 = \sum_n C_n \lambda_n X_n$ 

$Y_{3}$ =  $A Y_2 = \sum_n C_n \lambda_n^2 X_n$ 

In general: $Y_{i+1}$ =  $A Y_i = \sum_n C_n \lambda_n^i X_n$ 

After sufficiently many applications of $A$, the largest absolute eigenvalue will dominate completely 

$\dfrac{Y_i^\dagger A Y_i}{Y_i^\dagger Y_i} \Rightarrow max |\lambda_n| $

for normalised vectors, $Y_i^\dagger A Y_i \Rightarrow max |\lambda_n| $

In contrast to the power iteration, the `inverse power iteration` calculates the lowest eigenvalue of matrix $A$ as: 

$$A X_n = \lambda_n X_n  \Rightarrow \lambda^{-1}X_n = A^{-1} X_n$$

where, power iteration of $A^{-1}$ gives the lowest eigenvalue. However, instead of determining the $A^{-1}$ the following matrix equation is solved 

$A Y_{i+1} = Y_i$

Here, one starts from the guess vector $Y_i$ and solves the matrix equation to get to the next vector $Y_{i+1}$, which further goes through the same procedure. After numerous iterations 

$Y_i^\dagger A^{-1} Y_i \Rightarrow \dfrac{1}{min |\lambda_n|} $; for normalised vectors 

The `Inverse Power Iteration with shift` is utilised to obtain the smallest eigenvalue to a shifted matrix $A$, such that 

$A \Rightarrow A - \eta I$, where $\eta$ is the scalar constant close to the nearest eigenvalue.

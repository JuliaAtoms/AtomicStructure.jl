# Density matrices

The density matrix operator is defined as

$$\begin{equation}
\hat{\rho}_{ij} \defd \ketbra{i}{j}
\end{equation}$$

If the wavefunction is approximated using [Slater
determinants](https://en.wikipedia.org/wiki/Slater_determinant), where
one Slater determinant of $N$ electrons is defined as

$$\begin{equation}
\Phi(1,2,3,...,N) =
\frac{1}{\sqrt{N!}}
\left|\begin{matrix}
\chi_1(1)&
\chi_2(2)&
...&
\chi_N(N)
\end{matrix}\right|
\end{equation}$$

where $1,2,3...,N$ denote the $N$ different electronic coordinates
(spatial and spin) and $\chi_i$ are the $N$ different spin-orbitals,
the _transition density matrix_ for all $N$ electrons, between two
Slater determinants $\Phi_A$ and $\Phi_B$ is given by

$$\begin{equation}
\begin{aligned}
\rho_N^{AB}(1,...N;1',...,N')
&=
N!\Phi_A(1,...,N)\conj{\Phi_B}(1',...,N')\\
&=
\left|\begin{matrix}
\rho_1(1,1')&...&\rho_1(1,N')\\
\vdots&\ddots&\vdots\\
\rho_1(N,1')&...&\rho_1(N,N')
\end{matrix}\right|.
\end{aligned}
\end{equation}$$

## References

- Per-Olov Löwdin (1955). Quantum Theory of Many-Particle
  Systems. I. Physical Interpretations by Means of Density Matrices,
  Natural Spin-Orbitals, and Convergence Problems in the Method of
  Configurational Interaction. Physical Review, 97(6),
  1474–1489. [10.1103/physrev.97.1474](http://dx.doi.org/10.1103/physrev.97.1474)

- Per-Olov Löwdin (1955). Quantum Theory of Many-Particle
  Systems. II. Study of the Ordinary Hartree-Fock
  Approximation. Physical Review, 97(6),
  1490–1508. [10.1103/physrev.97.1490](http://dx.doi.org/10.1103/physrev.97.1490)

- Per-Olov Löwdin (1955). Quantum Theory of Many-Particle
  Systems. III. Extension of the Hartree-Fock Scheme to Include
  Degenerate Systems and Correlation Effects. Physical Review, 97(6),
  1509–1520. [10.1103/physrev.97.1509](http://dx.doi.org/10.1103/physrev.97.1509)


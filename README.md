# Nonvanishing a2 Coefficient

In the paper "Nonvanishing of the Second Coefficient of General Hecke Polynomials", we showed nonvanishing of the a2 coefficient for sufficiently large N,k.
This repository contains the code to compute the a2 coefficient in the finitely many remaining cases. 


- `compute_theta.py` is code to verify the claimed theta bounds from the paper. 
- `compute_a2_signs.sage` is code to compute a2(2,N,k). In particular, for any given m, it computes all pairs (N,k) for which a2(m,N,k) is positive/negative/zero.
- `final_a2_results_m2.txt` gives all pairs (N,k) for which a2(2,N,k) is positive or zero, along with dim S_k(Gamma_0(N)) and the actual value of a2(2,N,k). 
- `final_a2_results_m3.txt` gives all pairs (N,k) for which a2(3,N,k) is positive or zero, along with dim S_k(Gamma_0(N)) and the actual value of a2(3,N,k). 
- `final_a2_results_m4.txt` gives all pairs (N,k) for which a2(4,N,k) is negative or zero, along with dim S_k(Gamma_0(N)) and the actual value of a2(4,N,k). 

# SH_FWI
This is a  frequency domain 2-D SH wave  Full-Waveform Inversion (FWI) package implemented in OCTAVE. 9-Point finite difference stencil is used for forward solution. The implementation of the forward solution differs from that of another repository in how the coefficient matrix is built.The coefficient matrix is created without "nested for loops", which increase computation time in OCTAVE /MATLAB . First, all elements of this matrix are fully vectorized  and then combined into a sparse matrix.For inverse solution, the Limited-Memory BFGS approach is employed. A detailed description of the numerical methods used  will be uploaded.

(Note: The Pml_Damp_Chen function was created by using the Germania package, originally written in C by Daniel Köhn. In this OCTAVE version, "for loops" are not necessary. The misfit and gradient calculations are based on the SimpleFWI developed by Tristan van Leeuwen; however, unlike in SimpleFWI, the derivative matrix is used in the gradient calculation in this package. The optLBFGS-master was used to implement two loop recursion part of the Limited-Memory BFGS approach.)

# FORWARD SOLUTION
The article used for forwad solution  is  Chen, Z., Cheng, D.,Feng,W.,and Wu,T. 2013.An optimal 9-point finite difference  scheme for the Helmholtz equation with PML. International Journal Of Numerical  Analysis and Modeling,10,389-410.

# Perfectly matched layer Cheng et. al. (2013)

![chen](https://github.com/user-attachments/assets/a9194e57-c375-460d-aa18-aac144ccc8d2)


# Perfectly matched layer Lie et. al. (2023)

![Lie](https://github.com/user-attachments/assets/ea0a8560-52cc-47ca-8e5d-d7a88f4483b7)


# INVERSE SOLUTION

The main article used for inverse solution is Pratt, G.,R., Shin,C. and Hicks,G.,j.1998. Gauss-Newton and full Newton methods in frequency-space seismic waveform inversion. Geophysical Journal 

# True_S_Wave_Velocity_Model
<img width="753" height="381" alt="True Model" src="https://github.com/user-attachments/assets/b013d866-d75a-421c-962f-c249d128e51e" />

# Estimated_S_Wave_Velocity_Model

<img width="773" height="396" alt="Estimated_Model" src="https://github.com/user-attachments/assets/2d6a7859-f797-4b72-9395-8f8d731a21cf" />

# SH_FWI
This is a  frequency domain 2-D SH wave  Full-Waveform Inversion (FWI) package implemented in OCTAVE. 9-Point finite difference stencil is used for forward solution. The implementation of the forward solution differs from that of another repository in how the coefficient matrix is built.The coefficient matrix is created without "nested for loops", which increase computation time in OCTAVE /MATLAB . First, all elements of this matrix are fully vectorized  and then combined into a sparse matrix.For inverse solution, the Limited-Memory BFGS approach is employed. A detailed description of the numerical methods used  will be uploaded.

(Note: The Pml_Damp_Chen function was created by using the Germania package, originally written in C by Daniel Köhn. In this OCTAVE version, "for loops" are not necessary. The misfit and gradient calculations are based on the SimpleFWI developed by Tristan van Leeuwen; however, unlike in SimpleFWI, the derivative matrix is used in the gradient calculation in this package. The optLBFGS-master was used to implement two loop recursion part of the Limited-Memory BFGS approach.)

# FORWARD SOLUTION
# Nine Point Stencil
<img width="681" height="555" alt="image" src="https://github.com/user-attachments/assets/329896e3-d7ca-48ee-b8b3-e05ab22e0eb7" />

# 
  b ((A_4 (U_6-U_5 )-A_3 (U_5-U_4 ))/〖∆h〗^2 )                                                                                        +  (1-b)/2 ((A_2 (U_3-U_2 )-A_1 (U_2-U_1 ))/〖∆h〗^2 +(A_6 (U_9-U_8 )-A_5 (U_8-U_7 ))/〖∆h〗^2 )                      +b ((B_4 (U_8-U_5 )-B_3 (U_5-U_2 ))/〖∆h〗^2 )                                                                           +(1-b)/2 ((B_2 (U_7-U_4 )-B_1 (U_4-U_1 ))/〖∆h〗^2 +( B_6 (U_9-U_6 )-B_5 (U_6-U_3 ))/〖∆h〗^2 )+Cω^2 ρ(x,z)(cU_5+d/4 (U_8+U_6+U_2+U_4 )+e/4 (U_9+U_3+U_7+U_1 ))=f_y/〖∆h〗^2                                                                                                                                              

# Perfectly matched layer Cheng et. al. (2013)

![chen](https://github.com/user-attachments/assets/a9194e57-c375-460d-aa18-aac144ccc8d2)


# Perfectly matched layer Lie et. al. (2023)

![Lie](https://github.com/user-attachments/assets/ea0a8560-52cc-47ca-8e5d-d7a88f4483b7)


# INVERSE SOLUTION
# True_S_Wave_Velocity_Model
<img width="753" height="381" alt="True Model" src="https://github.com/user-attachments/assets/b013d866-d75a-421c-962f-c249d128e51e" />

# Estimated_S_Wave_Velocity_Model

<img width="773" height="396" alt="Estimated_Model" src="https://github.com/user-attachments/assets/2d6a7859-f797-4b72-9395-8f8d731a21cf" />

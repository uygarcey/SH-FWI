# SH_FWI
This is a  frequency domain 2-D SH wave  Full-Waveform Inversion (FWI) package implemented in OCTAVE. 9-Point finite difference stencil is used for forward solution. The implementation of the forward solution differs from that of another repository in how the coefficient matrix is built.The coefficient matrix is created without "nested for loops", which increase computation time in OCTAVE /MATLAB . First, all elements of this matrix are fully vectorized  and then combined into a sparse matrix.For inverse solution, the Limited-Memory BFGS approach is employed. A detailed description of the numerical methods used  will be uploaded.

(Note: The Pml_Damp_Chen function was created by using the Germania package, originally written in C by Daniel Köhn. In this OCTAVE version, "for loops" are not necessary. The misfit and gradient calculations are based on the SimpleFWI developed by Tristan van Leeuwen; however, unlike in SimpleFWI, the derivative matrix is used in the gradient calculation in this package. The optLBFGS-master was used to implement two loop recursion part of the Limited-Memory BFGS approach.)

# FORWARD SOLUTION
The article used for forwad solution  is  Chen, Z., Cheng, D.,Feng,W.,and Wu,T. 2013.An optimal 9-point finite difference  scheme for the Helmholtz equation with PML. International Journal Of Numerical  Analysis and Modeling,10,389-410.
Just try to fallow formulas and  derive it yourself. Descriptions will be provided later.

# 
 ![surfreq_sh](https://github.com/user-attachments/assets/dac4a8a3-08af-4071-92bd-569b848f0c74)  

 # 

 ![ekran](https://github.com/user-attachments/assets/590ed9f7-088f-4184-829e-ff9e7b4822d4)  

 # 

 ![2](https://github.com/user-attachments/assets/d463f2a5-8f38-4ee7-8859-755e277d3145)    

#
![3](https://github.com/user-attachments/assets/2eb77008-905e-422c-a367-f13dbd426115)      

#
![5](https://github.com/user-attachments/assets/badaf931-9949-4542-a7a4-5ea2a76b1bbe)

#
![6](https://github.com/user-attachments/assets/059b74bc-394d-4c2d-b16a-1c618723f2be)
#
![8](https://github.com/user-attachments/assets/52866555-2508-4e3c-a25a-f3d9ee884074)
#
![10](https://github.com/user-attachments/assets/2ea2507e-2904-4ed9-8677-491b1b344ee5)
#
![11](https://github.com/user-attachments/assets/932ca98b-3c20-4391-8eb3-74dc3b5b2f17)
#
![12](https://github.com/user-attachments/assets/3eb81755-6115-4631-914e-4b971f35d229)
#
![13](https://github.com/user-attachments/assets/778f00ab-4a43-494b-8def-e376ec16d374)
#
![14](https://github.com/user-attachments/assets/f175002f-81d7-46ab-9212-e525a0cc978b)
#
![15](https://github.com/user-attachments/assets/6981763e-ab8c-4f6b-bf57-e421207b90a5)
#

# Nine Point Stencil
<img width="681" height="555" alt="image" src="https://github.com/user-attachments/assets/329896e3-d7ca-48ee-b8b3-e05ab22e0eb7" />

# 
![16](https://github.com/user-attachments/assets/ffb23986-e978-4ee6-8754-888554343553)

#
The last thing is left. Derivation of cofeficenst corresponding each cell in 9-point stencil. It wll be soon..

                                                                                                                                          

# Perfectly matched layer Cheng et. al. (2013)

![chen](https://github.com/user-attachments/assets/a9194e57-c375-460d-aa18-aac144ccc8d2)


# Perfectly matched layer Lie et. al. (2023)

![Lie](https://github.com/user-attachments/assets/ea0a8560-52cc-47ca-8e5d-d7a88f4483b7)


# INVERSE SOLUTION
# True_S_Wave_Velocity_Model
<img width="753" height="381" alt="True Model" src="https://github.com/user-attachments/assets/b013d866-d75a-421c-962f-c249d128e51e" />

# Estimated_S_Wave_Velocity_Model

<img width="773" height="396" alt="Estimated_Model" src="https://github.com/user-attachments/assets/2d6a7859-f797-4b72-9395-8f8d731a21cf" />

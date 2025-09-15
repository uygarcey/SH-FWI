# SH_FWI
This is a  frequency domain 2-D SH wave  Full-Waveform Inversion (FWI) package implemented in OCTAVE. 9-Point finite difference stencil is used for forward solution. The implementation of the forward solution differs from that of another repository in how the coefficient matrix is built.The coefficient matrix is created without "nested for loops", which increase computation time in OCTAVE /MATLAB . First, all elements of this matrix are fully vectorized  and then combined into a sparse matrix.For inverse solution, the Limited-Memory BFGS approach is employed. A detailed description of the numerical methods used  will be uploaded.

(Note: The Pml_Damp_Chen function was created by using the Germania package, originally written in C by Daniel Köhn. In this OCTAVE version, "for loops" are not necessary. The misfit and gradient calculations are based on the SimpleFWI developed by Tristan van Leeuwen; however, unlike in SimpleFWI, the derivative matrix is used in the gradient calculation in this package. The optLBFGS-master was used to implement two loop recursion part of the Limited-Memory BFGS approach.)

# FORWARD SOLUTION
The article used for forwad solution  is  Chen, Z., Cheng, D.,Feng,W.,and Wu,T. 2013.An optimal 9-point finite difference  scheme for the Helmholtz equation with PML. International Journal Of Numerical  Analysis and Modeling,10,389-410.
Just try to fallow main formulas and  derive it yourself. Descriptions will be provided later.

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
The last thing is left. Derivation of cofeficenst corresponding each cell in 9-point stencil.

                                                                                                                                          

# Perfectly matched layer Cheng et. al. (2013)

![chen](https://github.com/user-attachments/assets/a9194e57-c375-460d-aa18-aac144ccc8d2)


# Perfectly matched layer Lie et. al. (2023)

![Lie](https://github.com/user-attachments/assets/ea0a8560-52cc-47ca-8e5d-d7a88f4483b7)


# INVERSE SOLUTION

The main article used for inverse solution is Pratt, G.,R., Shin,C. and Hicks,G.,j.1998. Gauss-Newton and full Newton methods in frequency-space seismic waveform inversion. Geophysical Journal International,133(2),341-362.
Just try to fallow main formulas and  derive it yourself. Descriptions will be provided later.
#
![İnverse_1](https://github.com/user-attachments/assets/8cb38c71-741d-4477-9576-2a33c6ace53f)

#
![Inverse_2](https://github.com/user-attachments/assets/237abc54-4bab-4ffa-93cb-ab59117dcd9d)

#
![Inverse_3](https://github.com/user-attachments/assets/53eba224-4ec0-4951-81e9-3897d78477fb)
#
![test](https://github.com/user-attachments/assets/d91628e7-f931-4118-8b6b-6497ddbd184e)


#
![Inverse_4](https://github.com/user-attachments/assets/84202507-a324-42a0-b025-5ec87a1ca2a9)

#
![inverse_5](https://github.com/user-attachments/assets/bbe1519f-e199-4e6d-ad04-b4b6611218d9)

#
![inverse_6](https://github.com/user-attachments/assets/1fe0ec7c-ac08-4d2d-a38a-a796b318b2f9)

#
![inverse_7](https://github.com/user-attachments/assets/a0763a25-2cc6-4303-8595-ff805dba7312)

#
![inverse_8](https://github.com/user-attachments/assets/e2804778-439c-440b-89eb-c8657d970a1b)

#
![inverse_9](https://github.com/user-attachments/assets/809a3e12-9f94-4809-a692-83592e98bcdf)

#
![inverse_10](https://github.com/user-attachments/assets/785f01b5-11ea-4edc-b56b-c7ef43d1a9c1)

#
![inverse_11](https://github.com/user-attachments/assets/7e910aa7-a954-445c-b7d9-016d394cc960)

#
![inverse_12](https://github.com/user-attachments/assets/03da5841-9c8c-4cf2-bab9-544065fdc142)


#
![test_test](https://github.com/user-attachments/assets/9d05f586-fae1-4c19-a130-eb6e668c8fa3)


# L-BFGS Nocedal 1980
#
![H_1](https://github.com/user-attachments/assets/59d5bfe2-ad8e-44f2-949c-c8a0fd2f4cad)

#
![H_2](https://github.com/user-attachments/assets/3d6a6ac6-c333-459a-854e-284df765f7b1)

#
![H_3](https://github.com/user-attachments/assets/05952f8a-87bf-411c-98c4-e2143bfad26f)

#
![H_4](https://github.com/user-attachments/assets/d33d088d-9076-4098-bcfb-acf043de898e)

#
![H_5](https://github.com/user-attachments/assets/43c16058-013c-4d47-a808-f99ea1aea5ea)

#
![H_6](https://github.com/user-attachments/assets/3cd82240-b583-4b14-b89a-d4a57890c3e8)

#
![H_7](https://github.com/user-attachments/assets/6f3a2b82-f25e-499b-8ef5-2dd215594384)

#
![H_8](https://github.com/user-attachments/assets/d999821d-742f-4442-9b69-66ed8e16f55f)

#
![H_9](https://github.com/user-attachments/assets/049a8d64-246c-4750-b564-e885412a6167)

#
![H_10](https://github.com/user-attachments/assets/7a3b4664-2520-4e23-ab6a-8c9aebd97bb1)

#
![H_11](https://github.com/user-attachments/assets/bd094065-1c1b-4993-b58d-7510069f80fd)

#
![H_12](https://github.com/user-attachments/assets/b72f7472-b60a-46b0-b622-fbe1268a2b07)

#
![H_13](https://github.com/user-attachments/assets/e96716cc-c48b-43fa-b64b-82fa1f813636)

#
![H_14](https://github.com/user-attachments/assets/51d14729-a4aa-428b-b9e7-3d0870062448)

#
![H_15](https://github.com/user-attachments/assets/2128c2e1-79fa-461d-be78-01bc680bc6bd)

#
![H_16](https://github.com/user-attachments/assets/38a7305a-eb80-4821-9133-e992921a679f)

#
![H_17](https://github.com/user-attachments/assets/911f405a-077f-411b-96da-eb55bc7572e7)

#
![H_18](https://github.com/user-attachments/assets/5f059a48-a26f-4794-aa11-f3f37caa5ea7)

#
![H_19](https://github.com/user-attachments/assets/3cbe4aa7-fff1-4c2d-9545-f3f5360a024f)

#
![H_20](https://github.com/user-attachments/assets/5258936c-c723-4492-978e-ca7977214166)

#
![H_21](https://github.com/user-attachments/assets/52a6e70d-f383-4828-9584-350aae24729c)

#
![H_22](https://github.com/user-attachments/assets/3ef4492f-47c1-4283-96ad-c5a9a206804f)

#

![H_23](https://github.com/user-attachments/assets/04af081b-01fb-4e66-8111-cfe68061e74e)
#
![H_24](https://github.com/user-attachments/assets/3a5e8192-a996-4565-97f8-7d5e6553702a)
#
![H_25](https://github.com/user-attachments/assets/674d7684-9f63-4165-8d4e-d3a1e8eba05c)



#
The last things are wolfe conditions

# True_S_Wave_Velocity_Model
<img width="753" height="381" alt="True Model" src="https://github.com/user-attachments/assets/b013d866-d75a-421c-962f-c249d128e51e" />

# Estimated_S_Wave_Velocity_Model

<img width="773" height="396" alt="Estimated_Model" src="https://github.com/user-attachments/assets/2d6a7859-f797-4b72-9395-8f8d731a21cf" />

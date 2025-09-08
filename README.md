# SH_FWI
This is a  frequency domain 2-D SH wave  Full-Waveform Inversion (FWI) package implemented in OCTAVE. 9-Point stencil finite difference is used for forward solution.Implementation of forward solution differs from that of another repository in terms of building coefficient matrix.The coefficient matrix is created without "nested for loops", which increase computation time  in OCTAVE /MATLAB . First, all elements of this matrix is fully vectorized  and then are combined into a sparse matrix . For inverse solution, the limited-memory BFGS approach is employed . A detailed description of the numerical methods used  will be uploaded.

# True_S_Wave_Velocity_Model
<img width="753" height="381" alt="True Model" src="https://github.com/user-attachments/assets/b013d866-d75a-421c-962f-c249d128e51e" />R

# Estimated_S_Wave_Velocity_Model

<img width="773" height="396" alt="Estimated_Model" src="https://github.com/user-attachments/assets/2d6a7859-f797-4b72-9395-8f8d731a21cf" />

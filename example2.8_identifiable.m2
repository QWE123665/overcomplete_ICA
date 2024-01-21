restart
R=QQ[b_1..b_4,b_5,l_1,l_2,l_3,l_4,l_5,l_6]
A_1=matrix{toList{0,1,2,3}}
A_2=matrix{toList{3,9,14,1}}
A_3=matrix{toList{1,11,13,-89735/6339}}
A_4=matrix{toList{-27417/160871,282663/36181,17,19}}
A_5=matrix{toList{1,0,0,0}}
A_6=matrix{toList{0,1,0,0}}
A_7=matrix{toList{0,0,1,0}}
A_8=matrix{toList{0,0,0,1}}
for i from 1 to 8 do B_i= toList{A_i**A_i}
B=matrix{B_1,B_2,B_3,B_4,B_5,B_6,B_7,B_8}
-- B is A\odot A
C=transpose(B)
F=submatrix'(C,{0,4,5,8,9,10,12,13,14,15},{4,5,6,7})
v=matrix{{l_1},{l_2},{l_3},{l_4}}
u=toList{transpose(F*v)}
-- find points in Vcal_I\cap W by looking at the off-diagonal terms
I=ideal(u_0_0_0-b_1*b_2,u_0_1_0-b_1*b_3,u_0_2_0-b_1*b_4,u_0_3_0-b_2*b_3,u_0_4_0-b_2*b_4,u_0_5_0-b_3*b_4)
L=minimalPrimes I
n=# L
for i from 0 to (n-1) list eliminate({l_1,l_2,l_3,l_4},L_i)
-- all the solutions correspond to the existing columns of A, so the matrix is identifiable

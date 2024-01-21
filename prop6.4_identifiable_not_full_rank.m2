restart
R=QQ[b_1..b_4,b_5,l_1,l_2,l_3,l_4,l_5,l_6]
A_1=matrix{toList{0,1,2,3,0}}
A_2=matrix{toList{3,9,14,1,0}}
A_3=matrix{toList{1,11,13,-89735/6339,0}}
A_4=matrix{toList{-27417/160871,282663/36181,17,19,0}}
A_5=matrix{toList{1,1,1,1,1}}
A_6=matrix{toList{1,2,3,4,7}}
A_7=matrix{toList{1,0,0,0,0}}
A_8=matrix{toList{0,1,0,0,0}}
A_9=matrix{toList{0,0,1,0,0}}
A_10=matrix{toList{0,0,0,1,0}}
A_11=matrix{toList{0,0,0,0,1}}
for i from 1 to 11 do B_i= toList{A_i**A_i}


-- checking (A_1,I_5) identifable 
D_1=matrix{B_1,B_2,B_3,B_4,B_7,B_8,B_9,B_10,B_11}
-- B is A\odot A
C=transpose(D_1)
F=submatrix'(C,{0,5,6,10,11,12,15,16,17,18,20,21,22,23,24},{4,5,6,7,8})
v=matrix{{l_1},{l_2},{l_3},{l_4}}
u=toList{transpose(F*v)}
-- find points in Vcal_I\cap W by looking at the off-diagonal terms
I=ideal(u_0_0_0-b_1*b_2,u_0_1_0-b_1*b_3,u_0_2_0-b_1*b_4,u_0_3_0-b_1*b_5,u_0_4_0-b_2*b_3,u_0_5_0-b_2*b_4,u_0_6_0-b_2*b_5,u_0_7_0-b_3*b_4,u_0_8_0-b_3*b_5,u_0_9_0-b_4*b_5)
L=minimalPrimes I
n=# L
for i from 0 to (n-1) list eliminate({l_1,l_2,l_3,l_4},L_i)
-- all the solutions correspond to the existing columns of A, so the matrix is identifiable
rank D_1
-- rank is 8 so A\odot A does not have full rank


-- checking (A_1,a_5,I_5) identifable 
D_2=matrix{B_1,B_2,B_3,B_4,B_5,B_7,B_8,B_9,B_10,B_11}
-- B is A\odot A
C=transpose(B)
F=submatrix'(C,{0,5,6,10,11,12,15,16,17,18,20,21,22,23,24},{5,6,7,8,9})
v=matrix{{l_1},{l_2},{l_3},{l_4},{l_5}}
u=toList{transpose(F*v)}
-- find points in Vcal_I\cap W by looking at the off-diagonal terms
I=ideal(u_0_0_0-b_1*b_2,u_0_1_0-b_1*b_3,u_0_2_0-b_1*b_4,u_0_3_0-b_1*b_5,u_0_4_0-b_2*b_3,u_0_5_0-b_2*b_4,u_0_6_0-b_2*b_5,u_0_7_0-b_3*b_4,u_0_8_0-b_3*b_5,u_0_9_0-b_4*b_5)
L=minimalPrimes I
n=# L
for i from 0 to (n-1) list eliminate({l_1,l_2,l_3,l_4,l_5},L_i)
-- all the solutions correspond to the existing columns of A, so the matrix is identifiable
rank D_2
-- rank is 9 so A\odot A does not have full rank


-- checking (A_1,a_5,a_6,I_5) identifable 
D_3=matrix{B_1,B_2,B_3,B_4,B_5,B_6,B_7,B_8,B_9,B_10,B_11}
-- B is A\odot A
C=transpose(B)
F=submatrix'(C,{0,5,6,10,11,12,15,16,17,18,20,21,22,23,24},{6,7,8,9,10})
v=matrix{{l_1},{l_2},{l_3},{l_4},{l_5},{l_6}}
u=toList{transpose(F*v)}
-- find points in Vcal_I\cap W by looking at the off-diagonal terms
I=ideal(u_0_0_0-b_1*b_2,u_0_1_0-b_1*b_3,u_0_2_0-b_1*b_4,u_0_3_0-b_1*b_5,u_0_4_0-b_2*b_3,u_0_5_0-b_2*b_4,u_0_6_0-b_2*b_5,u_0_7_0-b_3*b_4,u_0_8_0-b_3*b_5,u_0_9_0-b_4*b_5)
L=minimalPrimes I
n=# L
for i from 0 to (n-1) list eliminate({l_1,l_2,l_3,l_4,l_5},L_i)
-- all the solutions correspond to the existing columns of A, so the matrix is identifiable
rank D_3
-- rank is 10 so A\odot A does not have full rank


-- checking rank A=4, J=7, then A\odot A does not have full rank unless it has collinear columns
--without loss of generality, A has size 4*7 and the last four columns of A form the identity matrix
restart
R=QQ[a_1..a_4,b_1..b_4,c_1..c_4,l_1..l_3]
A_1=matrix{toList{a_1..a_4}}
A_2=matrix{toList{b_1..b_4}}
A_3=matrix{toList{c_1..c_4}}
A_4=matrix{toList{1,0,0,0}}
A_5=matrix{toList{0,1,0,0}}
A_6=matrix{toList{0,0,1,0}}
A_7=matrix{toList{0,0,0,1}}
for i from 1 to 7 do B_i= toList{A_i**A_i}
B=matrix{B_1,B_2,B_3,B_4,B_5,B_6,B_7}
I=minors(7,B);
minimalPrimes I
-- all the minimal primes of I correspond to some collumns of A are collinear or there exist three columns of A forming a matrix of rank 2

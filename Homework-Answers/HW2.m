close all;
clear all;

A = [55 -10 -20
    -10 30 -15
    -20 -15 65];
[L,U,P] = lu(A);

concatenated_matrix = [A P L U];

save('A1.dat','concatenated_matrix','-ascii')

LU_Matrix = zeros(3,100);
LU_Matrix2 = zeros(3,100);
for j = 1:100
    b = [50;0;j];
    d = P*b;
    y = L\d;
    x = U\y;
    x_2 = inv(A)*b;
    LU_Matrix(:,j) = [x];
    LU_Matrix2(:,j) = [x_2];
end
LU_Matrix;
LU_Matrix2;
LU_Matrix3 = abs(LU_Matrix - LU_Matrix2);
save('A2.dat','LU_Matrix','-ascii')
save('A3.dat','LU_Matrix3','-ascii')

s = sqrt(2)/2;
A_2 = [-s,1,0,0,0,0,0,0,0,s,0,0,0
        -s,0,0,0,0,0,0,0,-1,-s,0,0,0
        0,-1,1,0,0,0,0,0,0,0,0,0,0
        0,0,0,0,0,0,0,0,0,0,-1,0,0
        0,0,-1,s,0,0,0,0,0,0,0,-s,0
        0,0,0,-s,0,0,0,0,0,0,0,-s,-1
        0,0,0,-s,-1,0,0,0,0,0,0,0,0
        0,0,0,0,1,-1,0,0,0,0,0,0,0
        0,0,0,0,0,0,0,0,0,0,0,0,1
        0,0,0,0,0,1,-1,0,0,-s,0,s,0
        0,0,0,0,0,0,0,0,0,s,1,s,0
        0,0,0,0,0,0,1,-1,0,0,0,0,0
        0,0,0,0,0,0,0,0,1,0,0,0,0];
b_2 = [0, 0, 0, 0, 0, 0, 0, 0, 5, 0, 5, 0, 5].';
[L,U,P] = lu(A_2);
permutation = P * b_2;
y_2 = L \ permutation;
x_2 = U \ y_2;
x_3 = A_2\b_2;
save('A4.dat','y_2','-ascii')
save('A5.dat','x_2','-ascii')
save('A6.dat','x_3','-ascii')
q = 0;
while abs(norm(x_2,Inf)) <= 30
    q = q + .01;
    b_2 = [0, 0, 0, 0, 0, 0, 0, 0, 5, 0, 5 + q, 0, 5].';
    x_2 = A_2\b_2;
end
weight = q+5;
save('A7.dat','weight','-ascii')

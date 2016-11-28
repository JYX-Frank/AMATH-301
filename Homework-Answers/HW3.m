clear all; close all; clc; format long;

A_50 = 2*diag(ones(50,1),0) - diag(ones(49,1),1) - diag(ones(49,1),-1);

ro = zeros(50,1);
for j = 1:50
    ro(j,:) = 2*(1-cos(23*pi/51))*sin((23*pi*j/51));
end

D = diag(diag(A_50));
R = A_50 - D;

% 1 A
A_1 = inv(D)*(D-A_50);
A_1(:,51) = inv(D)*ro;
save('A1.dat','A_1','-ascii');

guess = ones(50,1);

%Set tolerance
tol = 1e-4;
error = 2 * tol;
new_matrix = guess;
iteration = 1;
for i = 1:100000
    iteration = iteration + 1;
    v_new = D\(ro-R*guess);
    guess = v_new;
    new_matrix(:,i+1) = v_new;
    error = norm(new_matrix(:,iteration)-new_matrix(:,iteration-1),Inf);
    
    %compare and stop if tolerance is met
    if error <= tol
        break;
    end
  
end

trueIteration1 = iteration - 1;
endMat = new_matrix(:,end);
save('A3.dat','trueIteration1','-ascii');
save('A2.dat','endMat','-ascii');

S = tril(A_50);
T = triu(A_50) - D;

% 1 C
C_3 = inv(S) * - T;
C_3(:,51) = inv(S)*ro;

save('A4.dat','C_3','-ascii');
x0 = ones(50,1);
X(:,1) = x0;
iterations = 1;
error = 2*tol;
while error > tol
    iterations = iterations + 1;
    X(:,iterations) = S \ (ro-T*X(:,iterations-1));
    error = norm(A_50*X(:,iterations)-ro,Inf);
end

endMat2 = X(:,end);
save('A5.dat','endMat2','-ascii');
trueIterations = iterations - 1;
save('A6.dat','trueIterations','-ascii');
% 2A
L = tril(A_50) - D;
U = triu(A_50) - D;

M = inv(D+1.5*L) * (-(1.5*U+(1.5-1)*D));
M(:,51) = inv(D+1.5*L) * (1.5*ro);

save('A7.dat','M','-ascii');
%2B
listomegas = 1:0.01:1.99;
new_matrix2 = zeros(100,50);
for j = 1:100
    omega = listomegas(j);
    matrix_B = -(D+omega*L)\(omega*U+(omega-1)*D);
    new_matrix2(j,:) = abs(eig(matrix_B));
end
A_8 = new_matrix2;
save('A8.dat','A_8','-ascii');

smallestMatrix =  [min(max(new_matrix2.'))+1;
    min(max(new_matrix2.'))];

save('A9.dat','smallestMatrix','-ascii');

guessOfOnes = ones(50,1);
Y(:,1) = guessOfOnes;
w = 1:0.01:1.99;

for p = 1:100
    omega = w(p);
    phi = guessOfOnes;
    M_x = inv(D+(omega)*L) * (-((omega)*U+(omega-1)*D));
    C_x = inv(D+omega*L) * (omega*ro);
    
    for q = 1:200
        phi = M_x*phi + C_x;
    end
    
    resid(p) =norm(A_50*phi-ro,2);
end
resid2 = resid.';
save('A10.dat','resid2','-ascii');
tol = 1e-4;
error = 2 * tol;
it = 1;
phi_2 = guessOfOnes;
while error > tol
    it = it + 1;
    M_x = inv(D+(1.89)*L) * (-((1.89)*U+(1.89-1)*D));
    C_x = inv(D+1.89*L) * (1.89*ro);
    phi_2 = M_x * phi_2 + C_x;
    error = norm(A_50*phi_2-ro,Inf);
end
save('A11.dat','phi_2','-ascii');
save('A12.dat','it','-ascii');

clear all; close all; clc;
t = 1;
h = 0.1;
dt = 0:h:t;
%dx/dt = -10x + 1.5sin(pi*t)
x = zeros(1,length(dt));
x(1) = 0;
%1A
dx_dt = @(t,x) -2*(x/(.25+x))*x + 1.5*sin(pi*t);
[t_FW, x_FW] = Forward_Euler(dx_dt,1,1,0.1);
x_FW = x_FW.';
save('A1.dat','x_FW','-ascii');

%1B
[t_out1,x_out1] = ode45(@(t_0,x_0) -2*(x_0/(.25+x_0))*x_0+1.5*sin(pi*t_0),dt,1);
save('A2.dat','x_out1','-ascii');

%1c 
[t_FW2, x_FW2] = Forward_Euler(dx_dt,1,1,0.01);
x_FW2 = x_FW2.';
%1c

[t_out3,x_out3] = ode45(@(t_0,x_0) -2*(x_0/(.25+x_0))*x_0+1.5*sin(pi*t_0),0:0.01:1,1);

newX = [x_FW2,x_out3];
save('A3.dat','newX','-ascii');

%1D
[t_out2,x_out2] = ode45(@(t_0,x_0) -2*(x_0/(.25+x_0))*x_0+1.5*sin(pi*t_0),[0;1],1);

new_x = [t_out2,x_out2];
save('A4.dat','new_x','-ascii');

%2A
A = [0 1
    -10/8 -3];

save('A5.dat','A','-ascii');

%2B
val = 0.8;
save('A6.dat','val','-ascii');

%2c
[t_1,x_1] = linear_forward_euler(A,[1,0],50,val.^2);
x_1 = x_1(1,1:end);

save('A7.dat','x_1','-ascii');

%2D
[t_2,x_2] = linear_forward_euler(A,[1,0],50,1.05*val);

x_2 = x_2(1,1:end);
save('A8.dat','x_2','-ascii');

%2E

g = -10;
L = 8;
delta = 3;
x_doubleprime = @(t,x) A*x;
[t_e, A10] = ode45(x_doubleprime,[0,50],[1;0]);
x_e = A10(:,1);
v_e = A10(:,2);
save('A9.dat', 'x_e', '-ascii');
save('A10.dat', 'v_e', '-ascii');

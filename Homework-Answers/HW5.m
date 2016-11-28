clear all; close all;

theta = 0:0.05*pi:2*pi;
r = 0:0.5:20;
[TH R] = meshgrid(theta,r);
constant = (sqrt(2)/(81 * sqrt(pi)) );
F = (abs(constant.*(6.*R - R.^2) .* exp(-R/3) .* cos(TH))).^2;
[TH R] = pol2cart(TH,R);

save('A1.dat','TH','-ascii');
save('A2.dat','R','-ascii');
save('A3.dat','F','-ascii');

%B
f1 = fminsearch('f',[0 1]);
f2 = fminsearch('f',[0 10]);
f3 = fminsearch('f',[pi 1]);
f4 = fminsearch('f',[pi 10]);

save('A4.dat','f1','-ascii');
save('A5.dat','f2','-ascii');
save('A6.dat','f3','-ascii');
save('A7.dat','f4','-ascii');

A8 = zeros(2,2);
A9 = zeros(2,2);
A10 = zeros(2,2);
A11 = zeros(2,2);


A8 = DG_out([0,1]);
A9 = DG_out([0 10]);
A10 = DG_out([pi 1]);
A11 = DG_out([pi 10]);

A8G = zeros(2,1);
A9G = zeros(2,1);
A10G = zeros(2,1);
A11G = zeros(2,1);

A8G = G([0 1]);
A9G = G([0 10]);
A10G = G([pi 1]);
A11G = G([pi 10]);

ans8 = newton2(@(input) G(input),@(input) DG_out(input),[0 1], 1e-4);
ans9 = newton2(@(input) G(input),@(input) G(input),[0 10], 1e-4);

% z = 22.4S + 26.29M + 39.88L
f = [22.4
    26.29
    39.88];
A = [0.1 .23 .31
    .22 .25 .38
    -.1 -.23 -.31
    .37 .42 .65
    0 0 0.27
    -1 0 0
    0 -1 0
    0 0 -1];
b = [5
    7
    -4
    15
    2
    -2
    -3
    -1];
lin_sol = round(linprog(-f,A,b),0);
save('A12.dat','f','-ascii');
save('A13.dat','A','-ascii');
save('A14.dat','b','-ascii');
save('A15.dat','lin_sol','-ascii');

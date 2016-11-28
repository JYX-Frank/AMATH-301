clear all; close all;

t = -90:10:110;
N = [7.24; 9.64; 12.87; 17.07; 23.19; 31.44; 38.56; 50.19; 62.98;
    76.21; 92.23; 106.02; 123.20; 132.16; 151.33; 179.32; 203.30;
    226.54; 248.71; 281.42; 307.75];

length = length(N);

dfdx = zeros(length,1);

dx = t(2)-t(1);

dfdx(1) = (-3*N(1) + (4*N(2)) - N(3))/(2*dx);
for i = 2:length-1
    dfdx(i) = (N(i+1) - N(i-1)) /(2*dx);
end
dfdx(end) = (3*N(end) - (4*N(end-1)) + N(end-2))/(2*dx);
save('A1.dat','dfdx','-ascii');

r = [0.308;
    0.325;
    0.342
    0.359
    .376
    .393
    .410
    .427
    .444 
    .461
    .478];
T = [640
    794
    885
    943
    1034
    1064
    1114
    1152
    1204
    1222
    1239];
theta = 0.7051;

h = (0.478-0.308)/10;
area = 0;

for i = 1:5
    area = area + r(2*i - 1)*T(2*i-1) + 4*(r(2*i)*T(2*i)) + r(2*i + 1)*T(2*i + 1);
end
area = (h/3)*area*theta;

h = (0.478-0.308)/10;
area2 = 0;

for i = 1:5
    area2 = area2 + r(2*i - 1) + 4*r(2*i) + r(2*i + 1);
end

area2 = (h/3)*area2*theta;

T_ = area/area2;

save('A2.dat','area','-ascii');
save('A3.dat','area2','-ascii');
save('A4.dat','T_','-ascii');

trap1 = trapz(r.*T)*h*theta;
trap2 = trapz(r)*h*theta;
trapArea = trap1/trap2;

save('A5.dat','trap1','-ascii');
save('A6.dat','trap2','-ascii');
save('A7.dat','trapArea','-ascii');

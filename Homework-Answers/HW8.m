clear all; close all;
clc;

x = .1:0.05:.5;
y = 1.8:0.05:2.2;
z = 11.8:0.05:12.2;

[X Y Z] = meshgrid(x,y,z);

X = X(:);
Y = Y(:);
Z = Z(:);

A1 = [X.'
    Y.'
    Z.'];

save('A1.dat','A1','-ascii');

%Initialize constants
a = 16;
b = 5;
c = 10;
d = 6;
e = 18;
f = 0.05;
dt = 0.01;

%Create three dimensional array
t = 0:dt:7;
t = t(:);
Dim3 = zeros(length(t),3,729);

dy = @(t,y)[
    -a*y(1,:)+ (f*y(2,:).*y(3,:));
    c*y(2,:) - (d*y(1,:).*y(3,:));
    -b*y(3,:) + (e*y(2,:).^2);
];

for i = 1:701
    Dim3(i,:,:) = A1;
    A1 = rk4singlestep(dy,dt,t(i),A1);
end

A2 = squeeze(Dim3(:,1,:));
A3 = squeeze(Dim3(:,2,:));
A4 = squeeze(Dim3(:,3,:));

save('A2.dat','A2','-ascii');
save('A3.dat','A3','-ascii');
save('A4.dat','A4','-ascii');

A5 = zeros(31,4);
A5(1,:) = 0.7;
r = [2.5 3.2 3.52 4];
for i = 2:31
    A5(i,:) = r.*A5(i - 1,:).*(1-A5(i - 1,:));
end

save('A5.dat','A5','-ascii');

A6 = [1;2;4];
save('A6.dat','A6','-ascii');

r = 2.5:0.001:4;
A7 = zeros(501,1501);
A7(1,:) = 0.7;

for i = 2:501
    A7(i,:) = r.*A7(i-1,:).*(1-A7(i-1,:));
end

a7 = A7(402:end,:);
save('A7.dat','a7','-ascii');

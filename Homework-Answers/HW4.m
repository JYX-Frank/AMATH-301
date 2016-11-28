clear all; close all; clc;
load salmon_data.csv;
t=(1:length(salmon_data)).';
new_matrix = zeros(2,2);
Vector_R = zeros(2,1);

new_matrix(1,1) = sum(t.^2);
new_matrix(1,2) = sum(t);
new_matrix(2,1) = sum(t);
new_matrix(2,2) = 77;
Vector_R(1) = sum(t.*salmon_data);
Vector_R(2) = sum(salmon_data);

p = new_matrix \ Vector_R;
save('A1.dat','new_matrix','-ascii');
save('A2.dat','Vector_R','-ascii');
save('A3.dat','p','-ascii');
p_2 = polyfit(t,salmon_data,2);
p_5 = polyfit(t,salmon_data,5);
p_8 = polyfit(t,salmon_data,8);

save('A4.dat','p_2','-ascii');
save('A5.dat','p_5','-ascii');
save('A6.dat','p_8','-ascii');

pre_vector = zeros(3,1);
pre_vector(1) = polyval(p_2,78);
pre_vector(2) = polyval(p_5,78);
pre_vector(3) = polyval(p_8,78);

save('A7.dat','pre_vector','-ascii');
new_time = 1:4:77.';
coarse_salmon = salmon_data(new_time);
save('A8.dat','coarse_salmon','-ascii');

nearest_neighbor = interp1(new_time,coarse_salmon,t,'nearest');
linear_interp = interp1(new_time,coarse_salmon,t,'linear');
cubic_interp = interp1(new_time,coarse_salmon,t,'cubic');
spline = interp1(new_time,coarse_salmon,t,'spline');

save('A9.dat','nearest_neighbor','-ascii');
save('A10.dat','linear_interp','-ascii');
save('A11.dat','cubic_interp','-ascii');
save('A12.dat','spline','-ascii');

RMS = zeros(4,1);

RMS(1) = sqrt( (1/77) * sum((salmon_data-nearest_neighbor).^2) );
RMS(2) = sqrt( (1/77) * sum((salmon_data-linear_interp).^2) );
RMS(3) = sqrt( (1/77) * sum((salmon_data-cubic_interp).^2) );
RMS(4) = sqrt( (1/77) * sum((salmon_data-spline).^2) );
save('A13.dat','RMS','-ascii');

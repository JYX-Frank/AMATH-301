clear all;
close all;
clc;

%1
dt = 0.004;
dx = 0.01;
a = 0.1;
L = 1;

A = 1 - 2 * ((a^2 * dt)/dx^2);
B = (a^2*dt)/dx^2;
M = B*diag(ones(98,1),1) + B*diag(ones(98,1),-1) + A*diag(ones(99,1),0);

save('A1.dat','A','-ascii');
save('A2.dat','B','-ascii');
save('A3.dat','M','-ascii');

A4 = max(abs(eig(M)));

save('A4.dat','A4','-ascii');

phi = @(x) exp(-200.*(x-0.5).^2);
y = zeros(99,1);
y(1) = phi(0);

t = dx:dx:1-dx;

for i = 2:99
   y(i) = phi(t(i));
end

save('A5.dat','y','-ascii');

A6 = zeros(99,501);
A6(:,1) = y;
time = 0:dx:2;
for i = 1:500
    A6(:,i+1) = M*A6(:,i);
end

A6 = [zeros(1,501); A6; zeros(1,501)];

save('A6.dat','A6','-ascii');
%2
[V,fr] = audioread('noisy_message.wav');

fft_V = fft(V);

A7 = abs(fft_V(1:1000));
save('A7.dat','A7','-ascii');

%Vectorized
filtered = fft_V;
filtered(abs(filtered) < 50) = 0;

A8 = abs(filtered(1:1000));

save('A8.dat','A8','-ascii');

% Non-Vectorized
% for i = 1:length(fft_V)
%     if fft_V(i) < 50
%         fft_V(i) = 0;
%     end
% end

inv_fft = ifft(filtered);
%inv_fft = abs(inv_fft);
% sound(inv_fft,fr)

A9 = inv_fft(1:1000);
save('A9.dat','A9','-ascii');

new_V = zeros(108680,1);

for i = 1:13585:length(V)
    section = V(i:i+13584);
    fft_section = fft(section);
    fft_filtered = fft_section;
    fft_filtered(abs(fft_filtered) < 3.2*mean(abs(fft_section)) ) = 0;
    
    new_V(i:i+13584) = ifft(fft_filtered);
end

new_V = new_V(1:1000);

save('A10.dat','new_V','-ascii');

clear all
close all
clc

%% Set up
N = 1064;
x = randn(1, N);
noise_axis = 0:N-1;
b = 1; a = [1, 0.9]; %coefficients to filter the noise

temp = filter(b, a, x);
y = temp(41:end);
filtered_axis = 40:N-1;
clear temp

figure()
subplot(2,1,1)
plot(noise_axis, x)
title('Sequence x')
xlabel('samples (n)')
ylabel('Amplitude (Au)')
subplot(2,1,2)
plot(filtered_axis, y)
title('Sequence y')
xlabel('samples (n)')
ylabel('Amplitude (Au)')

%% Part 1
[h,w]=freqz(1,[1 0.9],512);
plot(w/(2*pi),abs(h).^2)
 
%% Part 2
addpath(genpath([fileparts(pwd), '\Subpart1']));

N = N-40;
estimated_per = pgm(y);
pgm_axis = 0:1/N:(N-1)/N;

figure()
plot(pgm_axis, estimated_per)
xlim([0 0.5])

%% Part 3
figure()
plot(pgm_axis, estimated_per)
xlim([0.4 0.5])




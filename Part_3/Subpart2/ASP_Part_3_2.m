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
figure(), hold on
plot(w/(2*pi),abs(h).^2)
% plot([-7 1], [-3 -3]) %to find the cut-off when changing to logarithmic scale
title('PSD estimate for AR(1) process')
xlabel('normalized frequency')
ylabel('Amplitude (Au)')
legend('Ideal PSD')

%% Part 2
addpath(genpath([fileparts(pwd), '\Subpart1']));

N = 1024;
estimated_per = pgm(y);
pgm_axis = 0:1/N:(N-1)/N;

figure(), hold on;
plot(w/(2*pi),abs(h).^2, 'LineWidth', 1.5)
plot(pgm_axis, estimated_per, 'r')
xlim([0 0.5])
title('PSD estimate for AR(1) process')
xlabel('normalized frequency')
ylabel('Amplitude (Au)')
legend('Ideal PSD', 'Estimated periodogram')
% axis([0 0.5 -inf inf])

%% Part 3
figure(), hold on;
plot(w/(2*pi),abs(h).^2, 'LineWidth', 1.5)
plot(pgm_axis, estimated_per, 'r')
xlim([0.4 0.5])
title('PSD estimate for AR(1) process')
xlabel('normalized frequency')
ylabel('Amplitude (Au)')
legend('Ideal PSD', 'Estimated periodogram')

%% Window
rect = zeros(1024, 1);
rect(462:512, 1) = 1;
sinc = fftshift(ifft(rect));
plot(abs(sinc));

sgangersber = conv(y, sinc);
pasgar = pgm(sgangersber);
plot(pasgar)

%% Part 4
ry = xcorr(y, 'unbiased');
acf_axis = -length(y)+1 : length(y)-1;
a1_est = -ry(1, ceil(end/2)+1)/ry(1, ceil(end/2));
var_est = ry(1, ceil(end/2)) + a1_est*ry(1, ceil(end/2)+1);

N = 1024;
freq = 0:1/N:(N-1)/N;
Py = var_est ./ abs(1 + a1_est * exp(-2*1i*pi.*freq)).^2;

[h_model,w_model]=freqz(var_est,[1 a1_est],512);
figure(), hold on
plot(pgm_axis, estimated_per, 'k') %generated with pgm
plot(w/(2*pi),abs(h).^2, 'b', 'LineWidth', 1.5)
plot(w_model/(2*pi),abs(h_model).^2, 'r', 'LineWidth', 1.5)
xlim([0 0.5])
title('PSD estimate for AR(1) process')
xlabel('normalized frequency')
ylabel('Amplitude (Au)')
legend('Estimated periodogram', 'Ideal PSD', 'Model based PSD')

% [h_model,w_model]=freqz(var_est,[1 a1_est],512);
% figure(), hold on
% plot(pgm_axis, estimated_per) %generated with pgm
% plot(w_model/(2*pi),abs(h_model).^2, 'r', 'LineWidth', 1.5)
% xlim([0 0.5])
% title('PSD estimate for AR(1) process')
% xlabel('normalized frequency')
% ylabel('Amplitude (Au)')
% legend('Estimated periodogram', 'Model based PSD')

figure()
plot(freq, Py) %hard coded estimated PSD from model

%% Bode 
num = 1;
den = [1 0.9];
trans = tf(num, den);
[mag, phase, wout] = bode(trans);

mag = squeeze(mag);                                             % Reduce (1x1xN) Matrix To (1xN)
phase= squeeze(phase);
magr2 = (mag/max(mag)).^2;                                      % Calculate Power Of Ratio Of ‘mag/max(mag)’
dB3 = interp1(magr2, [wout phase mag], 0.5, 'spline');          % Find Frequency & Phase & Amplitude of Half-Power (-3 dB) Point
figure(1)
subplot(2,1,1)
semilogx(wout, 20*log10(mag), '-b',  dB3(1), 20*log10(dB3(3)), '+r', 'MarkerSize',10)
grid
subplot(2,1,2)
semilogx(wout, phase, '-b',  dB3(1), dB3(2), '+r', 'MarkerSize',10)
grid



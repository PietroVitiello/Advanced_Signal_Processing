%% Part 5
clear all
close all
clc

load sunspot.dat
% sun = filter(1, [1 0.9], sunspot(:,2));
sun = sunspot(:,2);

normalize = 1;
if (normalize == 1)
    sun = detrend(sun);
end

% Ideal psd
[h_ideal,w_ideal]=freqz(1,[1 0.9],length(sun)/2);

%estimated psd using pgm
addpath(genpath([fileparts(pwd), '\Subpart1']));

N = length(sun);
estimated_per = pgm(sun);
pgm_axis = 0:1/N:(N-1)/N;

% model psd
ry = xcorr(sun, 'unbiased');
acf_axis = -length(sun)+1 : length(sun)-1;

for orders = 1:10
    [a_est, var_est] = aryule(sun, orders);
    [h(:, orders), w]=freqz(var_est^(1/2),a_est,N/2);
end

% ry = xcorr(sun, 'unbiased');
% acf_axis = -length(sun)+1 : length(sun)-1;
% a1_est = -ry(ceil(end/2)+1)/ry(ceil(end/2));
% var_est = ry(ceil(end/2)) + a1_est*ry(ceil(end/2)+1);

% freq = 0:1/N:(N-1)/N;
% Py = var_est ./ abs(1 + a1_est * exp(-2*1i*pi.*freq)).^2;

% [h,w]=freqz(var_est,[1 a1_est],512);

figure(), hold on;

plot(pgm_axis, estimated_per)
plot(w/(2*pi),abs(h(:,1)).^2)
plot(w/(2*pi),abs(h(:,2)).^2)
plot(w/(2*pi),abs(h(:,10)).^2)

xlim([0 0.5])
title('PSD estimate of normalized sunspot data')
ylabel('Amplitude (Au)')
xlabel('normalized frequency')
legend('Estimated periodogram', 'order 1 model based estimated PSD', 'order 2 model based estimated PSD', 'order 10 model based estimated PSD')

%% Plot the individual psd
% figure() %ideal
% plot(w_ideal/(2*pi),abs(h_ideal).^2)
% title('ideal')

figure() %order 2
which_order = 2;
plot(w/(2*pi),abs(h(:,which_order)).^2)
title('model')

figure() %order 5
which_order = 5;
plot(w/(2*pi),abs(h(:,which_order)).^2)
title('model')

figure() %order 10
which_order = 10;
plot(w/(2*pi),abs(h(:,which_order)).^2)
title('model')

figure()
plot(pgm_axis, estimated_per)
xlim([0 0.5])
title('pgm')

%% pgm of not filtered sunspot
sgar = pgm(sunspot(:,2));
figure()
plot (sgar)


%% Part 5
clear all
close all
clc

load sunspot.dat
sun = filter(1, [1 0.9], sunspot(:,2));

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
    [a_est(orders, 1:orders+1), var_est(orders)] = aryule(sun, orders);
    [h(:, orders), w(:,orders)]=freqz(var_est(orders),[a_est(orders,1:orders+1)],N/2);
end

% ry = xcorr(sun, 'unbiased');
% acf_axis = -length(sun)+1 : length(sun)-1;
% a1_est = -ry(ceil(end/2)+1)/ry(ceil(end/2));
% var_est = ry(ceil(end/2)) + a1_est*ry(ceil(end/2)+1);

% freq = 0:1/N:(N-1)/N;
% Py = var_est ./ abs(1 + a1_est * exp(-2*1i*pi.*freq)).^2;

% [h,w]=freqz(var_est,[1 a1_est],512);

% Plot the various psd
figure()
plot(w_ideal/(2*pi),abs(h_ideal).^2)
title('ideal')

figure()
which_order = 10;
plot(w(:,which_order)/(2*pi),abs(h(:,which_order)).^2)
title('model')

figure()
plot(pgm_axis, estimated_per)
title('pgm')

%% pgm of not filtered sunspot
sgar = pgm(sunspot(:,2));
figure()
plot (sgar)


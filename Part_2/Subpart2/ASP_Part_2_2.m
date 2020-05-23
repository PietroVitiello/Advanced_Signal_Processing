%% Noise set up
clear all, close all
x = randn(1, 1000);
y = filter(2*ones(9,1), [1], x);

[ccf, lags] = xcorr(x, y, 'unbiased');

%% Part 1 full ACF
figure()
plot(lags, ccf);
title('Cross-correlation of uncorrelated process and its filtered version')
ylabel('Amplitude (Au)')
xlabel('relative distance between samples ({\tau})')
axis([-999 999 -inf inf])
%stem(lags, acf);

%% Part 1 stem between -20 and 20
figure()
stem(lags, ccf);
title('Zoomed Cross-correlation of uncorrelated process and its filtered version')
ylabel('Amplitude (Au)')
xlabel('relative distance between samples ({\tau})')
axis([-20 20 -inf inf])
%stem(lags, acf);

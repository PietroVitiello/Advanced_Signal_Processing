%% Noise set up
clear all
wnoise = randn(1, 1000);

[acf, lags] = xcorr(wnoise, 'unbiased');

%% Part 1
plot(lags, acf);
title('ACF of 1000 samples white Gaussian noise')
ylabel('Amplitude (Au)')
xlabel('relative distance between samples ({\tau})')
axis([-999 999 -inf inf])
%stem(lags, acf);

%% Part 2
figure()
plot(lags, acf);
title('Zoomed ACF of 1000 samples white Gaussian noise')
ylabel('Amplitude (Au)')
xlabel('relative distance between samples ({\tau})')
axis([-50 50 -inf inf])

figure()
plot(lags, acf);
title('Zoomed ACF of 1000 samples white Gaussian noise')
ylabel('Amplitude (Au)')
xlabel('relative distance between samples ({\tau})')
axis([600 inf -inf inf])

%% Part 4
filtered_wnoise = filter(2*ones(9,1), [1], wnoise);

%% Part 4 plot filtered noise
figure()
subplot(2, 1, 1)
plot(1:1000, wnoise)
subplot(2, 1, 2)
plot(1:1000, filtered_wnoise)

%% Part 4 plot filtered noise ACF
[f_acf, f_lags] = xcorr(filtered_wnoise, 'unbiased');
figure()
stem(f_lags, f_acf);
title('Autocorrelation function of filtered white Gaussian noise')
ylabel('Amplitude (Au)')
xlabel('relative distance between samples ({\tau})')
axis([-20 20 -inf inf])
% xlim([-20 20])


%% Noise set up
wnoise = randn(1, 1000);

acf = xcorr(wnoise, 'unbiased');
tao = [-999, 999];

%% Part 1
%plot(lags, acf);
title('Autocorrelation function of white Gaussian noise')
ylabel('Amplitude (Au)')
xlabel('relative distance between samples ({\tau})')
%xlim([-999 999])
stem(lags, acf);

%% Part 2
plot(lags, acf);
title('Autocorrelation function of white Gaussian noise')
ylabel('Amplitude (Au)')
xlabel('relative distance between samples ({\tau})')
xlim([-50 50])

%% Part 4
filtered_wnoise = filter(ones(9,1), 1, wnoise);
figure()
stem(lags, acf);
title('Autocorrelation function of filtered white Gaussian noise')
ylabel('Amplitude (Au)')
xlabel('relative distance between samples ({\tau})')
xlim([-20 20])


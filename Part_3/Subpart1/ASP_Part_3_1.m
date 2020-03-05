%% Set up
clc
clear all
close all

N = 128;
test_noise = randn(1, N);
test_axis = 0:1/N:(N-1)/N;

test_periodogram = pgm(test_noise);
disp(mean(test_periodogram))

figure()
plot(test_axis, test_periodogram)

%% Part 1
clear all

N = 128;
noise = randn(1, N);
per_axis = 0:1/N:(N-1)/N;

per = pgm(noise);
filtered_per = filter(0.2*(ones(1, 5)), 1, per);

figure(), hold on
plot(per_axis, per)
plot(per_axis, filtered_per, '-r', 'MarkerSize', 1)
title('Smoothed Periodogram estimate')
xlabel('normalized frequecy')
ylabel('Amplitude (Au)')
legend('original periodogram', 'smoothed periodogram')

%% Part 2
clc

noise_length = 32768;
noise = randn(1, noise_length);
noise_windows = [];
per_windows = [];

width_window = 128;
assert(rem(noise_length, width_window) == 0, 'This window size cannot be applied to 1028 samples')

for i = 1:(noise_length/width_window)
    noise_windows(i, :) = noise(1, (i-1)*width_window+1 : (i-1)*width_window+width_window);
    per_windows(i, :) = pgm(noise_windows(i, :));
end

disp('Windows periodograms done')

%% Part 3
per_axis = 0:1/N:(N-1)/N;
averaged_per = mean(per_windows);
plot(per_axis, averaged_per)


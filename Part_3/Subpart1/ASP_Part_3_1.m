%% Set up
clc
clear all
close all

for N = [128 256 512 9000]
    test_noise = randn(1, N);
    test_axis = 0:1/N:(N-1)/N;

    test_periodogram = pgm(test_noise);
    fprintf('mean: %f\n', mean(test_periodogram))
    fprintf('variance: %f\n\n', var(test_periodogram))

    figure()
    plot(test_axis, test_periodogram)
    title('Periodogram estimate')
    xlabel('normalized frequecy')
    ylabel('Amplitude (Au)')
end

%% Part 1
clear all
close all
clc

for N = [128 256 512]
    noise = randn(1, N);
    per_axis = 0:1/N:(N-1)/N;

    per = pgm(noise);
    filtered_per = filter(0.2*(ones(1, 5)), 1, per);
    fprintf('smoothed mean: %f\n', mean(filtered_per))
    fprintf('original mean: %f\n\n', mean(per))
    fprintf('variance: %f\n', var(filtered_per))
    fprintf('original variance: %f\n\n\n\n', var(per))

    figure(), hold on
    plot(per_axis, per)
    plot(per_axis, filtered_per, '-r', 'LineWidth', 1.5)
    title('Smoothed Periodogram estimate')
    xlabel('normalized frequecy')
    ylabel('Amplitude (Au)')
    legend('original periodogram', 'smoothed periodogram')
end

%% Part 2
clc

do_plot = 0;
titles = {'one' 'two' 'three' 'four' 'five' 'six' 'seven' 'eight'};

noise_length = 1024;
noise = randn(1, noise_length);
noise_windows = [];
per_windows = [];

width_window = 128;
per_axis = 0:1/width_window:(width_window-1)/width_window;
assert(rem(noise_length, width_window) == 0, 'This window size cannot be applied to 1028 samples')

for i = 1:(noise_length/width_window)
    noise_windows(i, :) = noise(1, (i-1)*width_window+1 : (i-1)*width_window+width_window);
    per_windows(i, :) = pgm(noise_windows(i, :));
    seg_variance(i) = var(per_windows(i,:));
    if do_plot
        subplot(2, 4, i)
        plot(per_axis, per_windows(i, :))
        title(sprintf('segment %s', titles{i}))
        xlabel('normalized frequecy')
        ylabel('Amplitude (Au)')
    end
end

disp('Windows periodograms done')
fprintf('\n\nmean variance of each periodogram: %f\n', mean(seg_variance))

%% Part 3
close all

averaged_per = mean(per_windows);
fprintf('variance of the averaged periodogram: %f\n', var(averaged_per))
fprintf('ratio between the variances: %f\n', mean(seg_variance)/var(averaged_per))

figure()
stem(per_axis, averaged_per)
title('Periodogram estimate')
xlabel('normalized frequecy')
ylabel('Amplitude (Au)')


clear all
close all
clc

%% Set up
% Simulation parameters
N = 1000;
if_normalize = 0;
std_noise = 0.1;

b = [1 2 3 2 1];
a = [1];
order = length(b);

% Signal preparation
x = randn(1, N);
y = filter(b, a, x);
if if_normalize == 1
    y = zscore(y);
end

noise = std_noise*randn(1, N);
z = y + noise;
SNR = 10*log10(std(y)^2/std(noise)^2);

%%% Alternatively load the signals that I used for the report
% load signals.mat

%% Finding optimal parameters
do_figures = 0;
optimal_sc = 15;

for repetitions = [1 2]
    %%% Finding best threshold factor
    thresholds = linspace(1, 30, 40);
    for param = thresholds
        [~, error_vector, ~] = lms_gear(x', z, 0.4, optimal_sc, param, order);

        avg_sqr_error(find(thresholds == param)) = mean(error_vector.^2);
    end
    
    if do_figures == 1
        figure()
        plot(thresholds, avg_sqr_error)
        title('Squared estimate error as a function of the threshold')
        xlabel('Threshold')
        ylabel('Squared estimate error')
    end

    %%% Finding best scaling factor
    optimal_th = thresholds(find(avg_sqr_error == min(avg_sqr_error)));
    scalings = linspace(1, 40, 40);
    for param = scalings
        [~, error_vector, ~] = lms_gear(x', z, 0.4, param, optimal_th, order);

        avg_sqr_error(find(scalings == param)) = mean(error_vector.^2);
    end
    
    if do_figures == 1
        figure()
        plot(scalings, avg_sqr_error)
        title('Squared estimate error as a function of the scaling')
        xlabel('Scaling')
        ylabel('Squared estimate error')
    end

    optimal_sc = scalings(find(avg_sqr_error == min(avg_sqr_error)));
end

%% Optimisation
step_size = 0.4;

[y_estimate, error_vector, weights] = lms_gear(x', z, step_size, optimal_sc, optimal_th, order);

figure()
plot(1:length(y_estimate), error_vector.^2)
title(['Squared estimate error using gear shifting LMS for \mu = ' num2str(step_size)])
xlabel('Time (n)')
ylabel('Squared estimate error')
xlim([0 400])

figure(), hold on
plot(1:length(y_estimate), weights)
scatter([1000 1000 1000], [1 2 3], 'b', 'filled')
title(['Evolution of the weights using gear shifting LMS for \mu = ' num2str(step_size)])
xlabel('Time (n)')
ylabel('Estimated Coefficients')
legend('w0', 'w1', 'w2', 'w3', 'w4', 'actual coefficients')

%% Compare to standard LMS
addpath(genpath([fileparts(pwd), '\Subpart2']));

step_size = 0.4;
[y_estimate, error_vector, weights] = lms(x', z, step_size, order);

figure(), hold on
plot(5:length(y_estimate)+4, weights)
scatter([1000 1000 1000], [1 2 3], 'b', 'filled')
title(['Evolution of the weights for \mu = ' num2str(step_size)])
xlabel('Time (n)')
ylabel('Estimated Coefficients')
legend('w0', 'w1', 'w2', 'w3', 'w4', 'actual coefficients')




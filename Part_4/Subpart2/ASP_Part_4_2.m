clear all
close all
clc

%% Set up
% Simulation parameters
N = 1000;
which_part = 1;
if_normalize = 0;

b = [1 2 3 2 1];
a = [1];
order = length(b);

% Signal preparation
x = randn(1, N);
y = filter(b, a, x);
if if_normalize == 1
    y = zscore(y);
end

if which_part == 1
    std_noise = 0.1;
    step_size = [0.002 0.08 0.1 0.5];
elseif which_part == 2
    noise_variance = [0.1 4.06 10];
    std_noise = noise_variance.^(1/2);
    step_size = 0.01;
else
    disp('Part can only be either 2 or 3') %should also throw error later because there is no noise_variance
end

for i = 1:length(std_noise)
    noise = std_noise(i)*randn(1, N);
    z(i, :) = y + noise;
    
    SNR(i,1) = 10*log10(std(y)^2/std(noise)^2);
end

%% Optimisation
for i = 1:length(std_noise)
    for gain = step_size
        which = find(step_size == gain);
        [y_estimate, error_vector, weights{which}] = lms(x', z(i,:), gain, order);
        
        if which_part == 1
            avg_error(which) = mean(abs(weights{which}(:, end)-b'));
            end_weights(:,which) = weights{which}(:, end);
        elseif which_part == 2
            avg_error(i) = mean(abs(weights{which}(:, end)-b'));
            end_weights(:,i) = weights{which}(:, end);
        end
        
        figure(), hold on
        plot(1:length(y_estimate), weights{which})
        scatter([1000 1000 1000], [1 2 3], 'b', 'filled')
        title(['Evolution of the weights for \mu = ' num2str(gain)])
        xlabel('Time (n)')
        ylabel('Estimated Coefficients')
        legend('w0', 'w1', 'w2', 'w3', 'w4', 'actual coefficients')
    end
end 

%% Squared error in estimated signal
if which_part == 1
    step_size = linspace(0.01, 0.27, 30);
    for gain = [0.01 0.1]
        [y_estimate, error_vector, weights] = lms(x', z, gain, order);
        
        avg_sqr_error(find(step_size == gain)) = mean(error_vector.^2);
        
        figure(), hold on
        plot(1:length(y_estimate), error_vector.^2)
        title(['Squared estimate error for \mu = ' num2str(gain)])
        xlabel('Time (n)')
        ylabel('Squared estimate error')
        xlim([0 400])
    end
end

figure(20)
plot(step_size, avg_sqr_error)
title('Squared estimate error as a function of the adaptation gain')
xlabel('Adaptation gain \mu')
ylabel('Squared estimate error')

%% There are too many variables
keep_vars = {'avg_error', 'SNR', 'weights', 'end_weights'};
clearvars ('-except', keep_vars{:})
    
    
    
    
clear all
close all
clc

%% Part 1
% Simulation parameters
N = 1000;
order = 5;
which_part = 1;
if_normalize = 0;

b = [1 2 3 2 1];
a = [1];

% Signal preparation
x = randn(1, N);
y = filter(b, a, x);
if if_normalize == 1
    y = zscore(y);
end

if which_part == 1
    std_noise = 0.1;
elseif which_part == 2
    noise_variance = linspace(0.1, 10, 6);
    std_noise = noise_variance.^(1/2);
else
    disp('Part can only be either 2 or 3') %should throw error ater because there is no noise_variance
end

for i = 1:length(std_noise)
    noise = std_noise(i)*randn(1, N);
    z(i, :) = y + noise;
    
    SNR(i,1) = 10*log10(std(y)^2/std(noise)^2);
end

step_size = [0.002 0.02 0.12 0.25 0.4 0.5];
for gain = step_size %linspace(0.002, 0.5, 6)
    [y_estimate, error_vector, weights{find(step_size == gain)}] = lms(x', z, gain, order);
    
    figure(), hold on
    plot(1:length(y_estimate), weights{find(step_size == gain)})
    scatter([1000 1000 1000], [1 2 3], 'b', 'filled')
    title(['Evolution of the weights for \mu = ' num2str(gain)])
    xlabel('Time (n)')
    ylabel('Estimated Coefficients')
    legend('w0', 'w1', 'w2', 'w3', 'w4')
end
    
    
    
    
    
clear all
close all
clc

%% Synthesis
% Simulation parameters
N = 1000;
if_normalize = 0;

a = [1 0.9 0.2];
b = [1];
order = length(a) - 1;

% Signal preparation
wgn = randn(1, N);
x = filter(b, a, wgn);
if if_normalize == 1
    x = zscore(x);
end

%% Analysis
gains = [0.002 0.01 0.07 0.2];

for step_size = gains
    [y_estimate, error, coeff] = adaptation_filter(x', step_size, order);

    figure(), hold on
    plot(1:length(y_estimate), coeff)
    scatter([1000 1000], -a(1, 2:3), 'b', 'filled')
    title(['Evolution of the coefficients for \mu = ' num2str(step_size)])
    xlabel('Time (n)')
    ylabel('Estimated Coefficients')
    legend('a1', 'a2', 'actual coefficients')
end





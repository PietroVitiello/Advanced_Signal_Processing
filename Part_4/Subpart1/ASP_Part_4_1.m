clear all
close all
clc

%% Set up

% Simulation parameters
N = 1000;
N_w = 4;
which_part = 2;
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
    
    SNR(i,1) = 20*log10(std(y)^2/std(noise)^2);
end

%% Weight optimisation
%finding R_xx
R_xx = zeros(N_w+1, N_w+1);
acf = xcorr(x, 'unbiased');

for i = 1:N_w+1
    R_xx(i, :) = acf(1, N+i-1 :-1: N+i-N_w-1);
end

for i = 1:length(std_noise)
    % finding p_zx
    p_zx = xcorr(z(i,:), x, 'unbiased');
    p_zx = p_zx(N : N+N_w)';
    
    %finding the weights for each different variance
    w_optimal(:, i) = inv(R_xx) * p_zx;
    error(i,1) = mean(abs(w_optimal(:, i) - b'));
end







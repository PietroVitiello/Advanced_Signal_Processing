close all
clear all
clc

addpath(genpath([fileparts(pwd), '\Subpart4']));
addpath(genpath([fileparts(pwd), '\Subpart5']));

%% Section 4.4
do_figures = 0;
repetitions = 400;

for i = 1:repetitions
    % Simulation parameters
    N = 6000;
    if_normalize = 0;
    step_size = 0.01;

    a = [1 0.9 0.2];
    b = [1];
    order = length(a) - 1;

    % Signal preparation
    wgn = randn(1, N);
    x = filter(b, a, wgn);
    if if_normalize == 1
        x = zscore(x);
    end

    %Otherwise it is possible to load the signal I used for the report
    % load signal.mat

    % Standard
    [y_estimate, ~, coeff] = adaptation_filter(x', step_size, order);
    variances(i, 1) = var(coeff(1500:end));

    if do_figures == 1
        figure(), hold on
        plot(1:length(y_estimate), coeff)
        scatter([N N], -a(1, 2:3), 'b', 'filled')
        title(['Evolution of the coefficients using Standard LMS'])
        xlabel('Time (n)')
        ylabel('Estimated Coefficients')
        legend('a1', 'a2', 'actual coefficients')
    end

    % Signed error
    [y_estimate, ~, coeff] = lms_signed_error(x', step_size, order);
    variances(i, 2) = var(coeff(1500:end));

    if do_figures == 1
        figure(), hold on
        plot(1:length(y_estimate), coeff)
        scatter([N N], -a(1, 2:3), 'b', 'filled')
        title(['Evolution of the coefficients using Signed error LMS'])
        xlabel('Time (n)')
        ylabel('Estimated Coefficients')
        legend('a1', 'a2', 'actual coefficients')
    end

    % Signed regressor
    [y_estimate, ~, coeff] = lms_signed_regressor(x', step_size, order);
    variances(i, 3) = var(coeff(1500:end));

    if do_figures == 1
        figure(), hold on
        plot(1:length(y_estimate), coeff)
        scatter([N N], -a(1, 2:3), 'b', 'filled')
        title(['Evolution of the coefficients using Signed regressor LMS'])
        xlabel('Time (n)')
        ylabel('Estimated Coefficients')
        legend('a1', 'a2', 'actual coefficients')
    end

    % Both signed
    [y_estimate, ~, coeff] = lms_sign_sign(x', step_size, order);
    variances(i, 4) = var(coeff(1500:end));

    if do_figures == 1
        figure(), hold on
        plot(1:length(y_estimate), coeff)
        scatter([N N], -a(1, 2:3), 'b', 'filled')
        title(['Evolution of the coefficients using Sign-Sign LMS'])
        xlabel('Time (n)')
        ylabel('Estimated Coefficients')
        legend('a1', 'a2', 'actual coefficients')
    end
end

average_variance = mean(variances);
clear vaiances

%% Section 4.5
fs = 44100;
N = 3000;
delay_samples = 1;
delay = zeros(delay_samples, 1);

[A, fs_rec] = audioread('A.m4a');
A = A(0.98*10^5 : 1.3*10^5, 1);
A = resample(A, fs, fs_rec);
boundary = 3400;
signal = A(boundary : boundary+N-1, 1);
delay_signal = [delay; signal(1: end-delay_samples)];

step_size = 0.7074;
order = 4;
do_figures = 2;

i = 3;
if i == 1
    e = zeros(1, 4);
end

x_low_lim = 0;
x_hig_lim = 1000;

% Standard
[y_estimate, ~, coeff] = adaptation_filter(delay_signal, step_size, order);
final_coeff(:, 1) = coeff(:, end);

if do_figures == 2
    pred = prediction(delay_signal, coeff(:, end), 0);
    e(i, 1) = mean((pred' - signal(order:end-1)).^2);
    figure(), hold on
    plot(1:length(pred), pred)
    plot(1:length(pred), signal(order:end-1))
    title(['Comparison between the prediction and the actual signal'])
    xlabel('Samples (n)')
    ylabel('Amplitude (Au)')
    legend('Prediction', 'Not delayed signal')
    xlim([x_low_lim x_hig_lim])
else
    figure()
    plot(1:length(y_estimate), coeff)
    title(['Evolution of the coefficients using Standard LMS'])
    xlabel('Time (n)')
    ylabel('Estimated Coefficients')
    legend('a1', 'a2', 'a3', 'a4')
end

% Signed error
[y_estimate, ~, coeff] = lms_signed_error(delay_signal, step_size, order);
final_coeff(:, 2) = coeff(:, end);

if do_figures == 2
    pred = prediction(delay_signal, coeff(:, end), 0);
    e(i, 2) = mean((pred' - signal(order:end-1)).^2);
    figure(), hold on
    plot(1:length(pred), pred)
    plot(1:length(pred), signal(order:end-1))
    title(['Comparison between the prediction and the actual signal'])
    xlabel('Samples (n)')
    ylabel('Amplitude (Au)')
    legend('Prediction', 'Not delayed signal')
    xlim([x_low_lim x_hig_lim])
else
    figure()
    plot(1:length(y_estimate), coeff)
    title(['Evolution of the coefficients using Signed error LMS'])
    xlabel('Time (n)')
    ylabel('Estimated Coefficients')
    legend('a1', 'a2', 'a3', 'a4')
end

% Signed regressor
[y_estimate, ~, coeff] = lms_signed_regressor(delay_signal, step_size, order);
final_coeff(:, 3) = coeff(:, end);

if do_figures == 2
    pred = prediction(delay_signal, coeff(:, end), 0);
    e(i, 3) = mean((pred' - signal(order:end-1)).^2);
    figure(), hold on
    plot(1:length(pred), pred)
    plot(1:length(pred), signal(order:end-1))
    title(['Comparison between the prediction and the actual signal'])
    xlabel('Samples (n)')
    ylabel('Amplitude (Au)')
    legend('Prediction', 'Not delayed signal')
    xlim([x_low_lim x_hig_lim])
else
    figure()
    plot(1:length(y_estimate), coeff)
    title(['Evolution of the coefficients using Signed regressor LMS'])
    xlabel('Time (n)')
    ylabel('Estimated Coefficients')
    legend('a1', 'a2', 'a3', 'a4')
end

% Both signed
[y_estimate, ~, coeff] = lms_sign_sign(delay_signal, step_size, order);
final_coeff(:, 4) = coeff(:, end);

if do_figures == 2
    pred = prediction(delay_signal, coeff(:, end), 0);
    e(i, 4) = mean((pred' - signal(order:end-1)).^2);
    figure(), hold on
    plot(1:length(pred), pred)
    plot(1:length(pred), signal(order:end-1))
    title(['Comparison between the prediction and the actual signal'])
    xlabel('Samples (n)')
    ylabel('Amplitude (Au)')
    legend('Prediction', 'Not delayed signal')
    xlim([x_low_lim x_hig_lim])
else
    figure()
    plot(1:length(y_estimate), coeff)
    title(['Evolution of the coefficients using Sign-Sign LMS'])
    xlabel('Time (n)')
    ylabel('Estimated Coefficients')
    legend('a1', 'a2', 'a3', 'a4')
end




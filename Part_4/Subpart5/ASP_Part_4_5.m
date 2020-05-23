clear all
close all
clc

addpath(genpath([fileparts(pwd), '\Subpart4']));

%% Setting up signals
fs = 44100;
N = 1000;
which_part = 1;
delay_samples = 1;
delay = zeros(delay_samples, 1);

[A,fs_rec] = audioread('A.m4a');
A = A(0.98*10^5 : 1.3*10^5, 1);
A = resample(A, fs, fs_rec);
boundary = 3400;
A = A(boundary : end, 1);

[E, ~] = audioread('E.m4a');
E = E(1.09*10^5 : 1.33*10^5, 1);
E = resample(E, fs, fs_rec);
boundary = 2700;
E = E(boundary : end, 1);

[S, ~] = audioread('S.m4a');
S = resample(S, fs, fs_rec);
boundary = 75440;
S = S(boundary : end, 1);

[T, ~] = audioread('T.m4a');
T = resample(T, fs, fs_rec);
boundary = 91300;
T = T(boundary : end, 1);

[X, ~] = audioread('X.m4a');
X = resample(X, fs, fs_rec);
boundary = 75560;
X = X(boundary : end, 1);

if which_part == 2
    new_fs = 16000;
    A = resample(A, new_fs, fs);
    E = resample(E, new_fs, fs);
    S = resample(S, new_fs, fs);
    T = resample(T, new_fs, fs);
    X = resample(X, new_fs, fs);
end

A = A(1 : N);
E = E(1 : N);
S = S(1 : N);
T = T(1 : N);
X = X(1 : N);

delay_A = [delay; A(1: end-delay_samples)];
delay_E = [delay; E(1: end-delay_samples)];
delay_S = [delay; S(1: end-delay_samples)];
delay_T = [delay; T(1: end-delay_samples)];
delay_X = [delay; X(1: end-delay_samples)];

letters = ['A', 'E', 'S', 'T', 'X'];
signal = {A, E, S, T, X};
delay_signal = {delay_A, delay_E, delay_S, delay_T, delay_X};

%% Finding the best parameters
best_order = 10*ones(5, 1);

for repetitions = [1 2]
    %%% Error for various gains
    tests = 100;
    steps = linspace(0.001, 1, tests);

    for ii = 1:length(letters)
        e = zeros(1, tests);

        for i = 1:tests
            [~, ~, coeff] = adaptation_filter(delay_signal{ii}, steps(i), best_order(ii));
            estimate = prediction(delay_signal{ii}, coeff(:, end), 0);
            e(1,i) = mean((estimate' - signal{ii}(best_order(ii):end-1)).^2);
        end

        [~, best_gain(ii, 1)] = min(e);
        best_gain(ii, 1) = steps(best_gain(ii, 1));
    end

    %%% Error for various orders
    p = 40;
    do_figures = 0;
    clear MDL AIC AIC_c

    for ii = 1:length(letters)
        e = zeros(1, p);

        for i = 1:p
            [~, ~, coeff] = adaptation_filter(delay_signal{ii}, best_gain(ii), i);
            estimate = prediction(delay_signal{ii}, coeff(:, end), 0);
            e(1,i) = sum((estimate' - signal{ii}(i:end-1)).^2);
        end

        MDL = log(e(1,:)) + (1:p).*log(N)./N;
        AIC = log(e(1,:)) + 2.*(1:p)./N;
        AIC_c = AIC + (2.*(1:p).*((1:p)+1))./(-(1:p)-1+N);
        cumulative = log(e(1,:));

        if do_figures == 1
            figure(), hold on;
            p_axis = 1:p;
            plot(p_axis, MDL)
            plot(p_axis, AIC)
            plot(p_axis, AIC_c)
            plot(p_axis, cumulative)
            legend('Minimum Description Length', 'Akaike Information Criterion', 'corrected AIC', 'Cumulative Error Squared')
            title(sprintf('Loss calculation for the letter "%s"', letters(ii)))
            xlabel("model order (p)")
            ylabel("Loss function")
        end

        [~, MDL_best(ii, 1)] = min(MDL);
        [~, AIC_best(ii, 1)] = min(AIC);
        [~, AICc_best(ii, 1)] = min(AIC_c);
        [~, cumulative_best(ii, 1)] = min(cumulative);
        
        best_order = AICc_best;

    end
end

%% Evolution of the coefficients
do_figures = 1;

for i = 1:length(letters)
    [signal_estimate, ~, coeff] = adaptation_filter(delay_signal{i}, best_gain(i), best_order(i));
    
    if do_figures == 1
        figure()
        plot(1:length(signal_estimate), coeff)
        title(['Evolution of the coefficients for the letter "', letters(i), '"'])
        xlabel('Time (n)')
        ylabel('Estimated Coefficient')
    end
    
    if do_figures == 2
        pred = prediction(delay_signal{i}, coeff(:, end), 0);
        figure(), hold on
        plot(1:length(pred), pred)
        plot(1:length(pred), signal{i}(best_order(i):end-1))
        title(['Comparison between the prediction and the actual signal for the letter "', letters(i), '"'])
        xlabel('Samples (n)')
        ylabel('Amplitude (Au)')
        legend('Prediction', 'Not delayed signal')
    end
end

%% Table of gain errors for Latex
gains = [0.2634 0.7074 0.8385 1];
e = zeros(4, 5);
for ii = 1:length(letters)
    for i = 1:length(sasa)
        [~, ~, coeff] = adaptation_filter(delay_signal{ii}, gains(i), best_order(ii));
        estimate = prediction(delay_signal{ii}, coeff(:, end), 0);
        e(i, ii) = mean((estimate' - signal{ii}(best_order(ii):end-1)).^2);
    end
end

%% Performance test
which_performance = 1;

%Comparing different orders
if which_performance == 1
    order_list = [1 4 10 29 33];
    performance = zeros(length(order_list), length(letters));

    for ii = 1:length(letters)
        for i = 1:length(order_list)
            [~, ~, coeff] = adaptation_filter(delay_signal{ii}, best_gain(ii), order_list(i));
            estimate = prediction(delay_signal{ii}, coeff(:, end), 0);
            e = estimate' - signal{ii}(order_list(i):end-1);

            performance(i, ii) = 10*log10(var(delay_signal{ii}) / var(e));
        end
    end
end

%for the optimal values
if which_performance == 2
    performance = zeros(1, length(letters));

    for ii = 1:length(letters)
            [~, ~, coeff] = adaptation_filter(delay_signal{ii}, best_gain(ii), best_order(ii));
            estimate = prediction(delay_signal{ii}, coeff(:, end), 0);
            e = estimate' - signal{ii}(best_order(ii):end-1);

            performance(1, ii) = 10*log10(var(signal{ii}(best_order(ii):end-1)) / var(e));
    end
end



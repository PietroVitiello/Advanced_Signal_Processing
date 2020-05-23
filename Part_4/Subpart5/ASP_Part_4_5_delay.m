clear all
close all
clc

addpath(genpath([fileparts(pwd), '\Subpart2']));
addpath(genpath([fileparts(pwd), '\Subpart4']));

%% Setting up signals
fs = 44100;
N = 1000;
delay_samples = 30;
delay = zeros(delay_samples, 1);

[A,fs_rec] = audioread('A.m4a');
A = A(0.98*10^5 : 1.3*10^5, 1);
A = resample(A, fs, fs_rec);
whole_a = A;
boundary = 3400;
A = A(boundary : boundary+N-1, 1);
delay_A = [delay; A(1: end-delay_samples)];

[E, ~] = audioread('E.m4a');
E = E(1.09*10^5 : 1.33*10^5, 1);
E = resample(E, fs, fs_rec);
boundary = 2700;
E = E(boundary : boundary+N-1, 1);
delay_E = [delay; E(1: end-delay_samples)];

[S, ~] = audioread('A.m4a');
S = S(0.98*10^5 : 1.32*10^5, 1);
S = resample(S, fs, fs_rec);
boundary = 3100;
S = S(boundary : boundary+N-1, 1);
delay_S = [delay; S(1: end-delay_samples)];

[T, ~] = audioread('A.m4a');
T = T(0.98*10^5 : 1.28*10^5, 1);
T = resample(T, fs, fs_rec);
boundary = 2450;
T = T(boundary : boundary+N-1, 1);
delay_T = [delay; T(1: end-delay_samples)];

[X, ~] = audioread('A.m4a');
X = X(0.98*10^5 : 1.27*10^5, 1);
X = resample(X, fs, fs_rec);
boundary = 4050;
X = X(boundary : boundary+N-1, 1);
delay_X = [delay; X(1: end-delay_samples)];

letters = ['A', 'E', 'S', 'T', 'X'];
signal = {A, E, S, T, X};
delay_signal = {delay_A, delay_E, delay_S, delay_T, delay_X};

%% Evolution of the coefficients
gain = 0.01;
order = 10;
do_figures = 0;

for i = 1:length(letters)
    [signal_estimate, e, coeff] = lms(delay_signal{i}, signal{i}, gain, order); %adaptation_filter(A, 0.01, 10);
    
    if do_figures == 1
        figure()
        plot(1:length(signal_estimate), coeff)
        title(['Evolution of the coefficients for the letter "', letters(i), '"'])
        xlabel('Time (n)')
        ylabel('Estimated Coefficient')
    end
end

figure(), hold on
nana = filter([coeff(:, end)'], 1, delay_A);
plot(nana)
% figure()
plot(A)

%% Error for various orders
p = 40;
do_figures = 1;
clear reference estimate MDL AIC AIC_c

for ii = 1:length(letters)
    e = zeros(1, p);
    reference = ones(N, p) .* signal{ii};
    
    for i = 1:p
        [~, ~, coeff] = lms(delay_signal{ii}, signal{ii}, gain, i);
        estimate(:,i) = filter([coeff(:, end)'], 1, delay_signal{ii});
%         estimate(:,i) = filter(1, [1 -coeff(:, end)'], delay_signal{ii});
    end

    for i = 1:N
        e(1,:) = (estimate(i,:) - reference(i,:)).^2 + e(1,:);
    end

    MDL = log(e(1,:)) + (1:p).*log(N)./N;
    AIC = log(e(1,:)) + 2.*(1:p)./N;
    AIC_c = AIC + (2.*(1:p).*((1:p)+1))./(-(1:p)-1+N);
    cumulative = log(e(1,:));

    %normalizing
    MDL = MDL/MDL(1);
    AIC = AIC/AIC(1);
    AIC_c = AIC_c/AIC_c(1);
    cumulative = cumulative/cumulative(1);
    
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
    [~, AICc_best(ii, 1)] = min(AIC_c);
    [~, cumulative_best(ii, 1)] = min(cumulative);
    
end

%% Performance test
order_list = [1 2 3 20 21 30 32];

for ii = 1:length(letters)
    for i = 1:length(order_list)
        [~, pa, coeff] = lms(delay_signal{ii}, signal{ii}, gain, order_list(i));
%         estimate = filter(1, [1 -coeff(:, end)'], delay_signal{ii});
        estimate = filter([coeff(:, end)'], 1, delay_signal{ii});
        
        e = estimate - signal{ii};
        performance (ii, i) = 10*log10(var(delay_signal{ii}) / var(e));
    end
end

% for ii = 1:length(letters)
%     for i = 1:length(order_list)
%         [~, e, ~] = lms(delay_signal{ii}, signal{ii}, gain, order_list(i));
%         
%         performance (ii, i) = 10*log10(var(delay_signal{ii}) / var(e));
%     end
% end

figure(), hold on
plot(estimate)
plot(delay_signal{ii})
% plot(signal{ii})



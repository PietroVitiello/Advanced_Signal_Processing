clc
clear all

load RRI_2_5.mat

RRI_fs = RRI_fs1;
trial = RRI_trial1;
clear RRI_fs1 RRI_fs3 RRI_fs3

figure()
plot(trial)
figure()
plot(RRI_trial2)
figure()
plot(RRI_trial3)

%get heart rate
h = 60 ./ trial;

%% Part 1
close all
clear he

alpha = [1 0.6];
he = [];

variance = var(h);
disp(variance)
disp(mean(h))

for n = [1 2]
    a = 1;
    for i = 0:10:(length(h)-8)
        if (i == 1000)
            he(a, n) = (1/8)* sum(alpha(n) .* h(i+1:end));
        else
            he(a, n) = (1/10)* sum(alpha(n) .* h(i+1:i+10));
        end
        a = a+1;
    end
end
variance = var(he);
disp(variance)
disp(mean(he))

n_bins = 20;

figure()
subplot(3,1,1)
histogram(h, n_bins, 'Normalization', 'pdf')
title('PDE of the original heart rate')
xlabel('BPM')
ylabel('Probability')
xlim([0 150])
subplot(3,1,2)
histogram(he(:,1), n_bins, 'Normalization', 'pdf')
title('PDE of the averaged heart rate for \alpha = 1')
xlabel('averaged BPM')
ylabel('Probability')
xlim([0 150])
subplot(3,1,3)
histogram(he(:,2), n_bins, 'Normalization', 'pdf')
title('PDE of the averaged heart rate for \alpha = 0.6')
xlabel('averaged BPM')
ylabel('Probability')
xlim([0 150])

%% Part 2
%close all

trial1 = detrend(RRI_trial1);
trial2 = detrend(RRI_trial2);
trial3 = detrend(RRI_trial3);

acf1 = xcorr(trial1, 'unbiased');
acf2 = xcorr(trial2, 'unbiased');
acf3 = xcorr(trial3, 'unbiased');

tau1 = -(length(trial1)-1) : length(trial1)-1;
tau2 = -(length(trial2)-1) : length(trial2)-1;
tau3 = -(length(trial3)-1) : length(trial3)-1;

figure()
subplot(3,1,1)
plot(tau1, acf1)
title('ACF of the first trial')
xlabel('lag')
ylabel('Amplitude (Au)')
subplot(3,1,2)
plot(tau2, acf2)
title('ACF of the second trial')
xlabel('lag')
ylabel('Amplitude (Au)')
subplot(3,1,3)
plot(tau3, acf3)
title('ACF of the third trial')
xlabel('lag')
ylabel('Amplitude (Au)')

% find order of the models
%% Partial Autocorrelation
close all

p = 10;

[a1, ~, k1] = aryule(trial1, p);
[a2, ~, k2] = aryule(trial2, p);
[a3, ~, k3] = aryule(trial3, p);

boundary1 = 1.96/sqrt(length(trial1));
boundary2 = 1.96/sqrt(length(trial2));
boundary3 = 1.96/sqrt(length(trial3));

orders = 1:p;

figure(1), hold on
stem(orders, -k1)
plot([1 p], [1 1]' * [boundary1 -boundary1], 'r--')
title('PAC function of the first trial')
xlabel('lag')
ylabel('Amplitude (Au)')

figure(2), hold on
stem(orders, -k2)
plot([1 p], [1 1]' * [boundary2 -boundary2], 'r--')
title('PAC function of the second trial')
xlabel('lag')
ylabel('Amplitude (Au)')

figure(3), hold on
stem(orders, -k3)
plot([1 p], [1 1]' * [boundary3 -boundary3], 'r--')
title('PAC function of the third trial')
xlabel('lag')
ylabel('Amplitude (Au)')

%% Loss Functions
clc

p = 10;
trials = {trial1' trial2' trial3'};
titles = {'first trial', 'second trial', 'third trial'};

for ii = 1:3
    clear filtered_trial a tri MDL AIC AIC_c
    N = length(trials{ii});
    tri = ones(N, p) .* trials{ii};

    filtered_trial = zeros(N, p);
    for i = 1:p
        a = aryule(trials{ii}, i);
        filtered_trial(:, i) = filter(-a, 1, trials{ii});
    end

    e = zeros(1, p);

    for i = 1:N
        e(1,:) = (filtered_trial(i,:) - tri(i,:)).^2 + e(1,:);
    end

    MDL = log(e(1,:)) + (1:p).*log(N)./N;
    AIC = log(e(1,:)) + 2.*(1:p)./N;
    AIC_c = AIC + (2.*(1:p).*((1:p)+1))./(-(1:p)-1+N);
%     cumulative = log(e(1,:));

    %normalizing
    MDL = MDL/MDL(1);
    AIC = AIC/AIC(1);
    AIC_c = AIC_c/AIC_c(1);
    
    if (ii == 2)
        MDL = -(MDL-1)+1;
        AIC = -(AIC-1)+1;
        AIC_c = -(AIC_c-1)+1;
    end
%     cumulative = cumulative/cumulative(1);

    figure(), hold on;
    p_axis = 1:p;
    plot(p_axis, MDL)
    plot(p_axis, AIC)
    plot(p_axis, AIC_c)
%     plot(p_axis, cumulative)
    legend('Minimum Description Length', 'Akaike Information Criterion', 'corrected AIC') %, 'Cumulative Error Squared')
    title(sprintf("Loss calculation for the %s", titles{ii}))
    xlabel("model order (p)")
    ylabel("Loss function")
end









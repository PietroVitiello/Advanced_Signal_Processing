clear all
close all
clc

%% Set up
do_normalization = 1;

load sunspot.dat
sun = sunspot(:,2);

if(do_normalization == 1)
    sun = zscore(sun);
end

%% Part 3
acf = xcorr(sun); % peak at 288
M = length(acf);

estimated_a = zeros(10);
for p = 1:10 % order
    if (p >= 1)
        padded_acf = [zeros(p,1); acf];
    else
        padded_acf = acf;
    end

    for i = 1:M
        for ii = 1:p
            H(i, ii) = padded_acf(p+i-ii,1);
        end
    end
    
    estimated_a(p,1:p) = inv(H'*H)*H'*acf;
end

estimated_a(2,1) = 1.414;
estimated_a(2,2) = -0.715;

%% Part 4
N = length(sun);
sun_matrix = ones(N, p) .* sun;

estimated_signal = zeros(N, p);
for i = 1:p
    estimated_signal(:, i) = filter([-1 estimated_a(i, :)], 1, sun);
end
figure(20), hold on
plot(estimated_signal(:,2))
plot(sun)

e = zeros(1, p);

for i = 1:N
    e(1,:) = (estimated_signal(i,:) - sun_matrix(i,:)).^2 + e(1,:);
end

MDL = log(e(1,:)) + (1:p).*log(N)./N;
AIC = log(e(1,:)) + 2.*(1:p)./N;
AIC_c = AIC + (2.*(1:p).*((1:p)+1))./(-(1:p)-1+N);
cumulative = log(e(1,:));

%normalizing
% MDL = MDL/MDL(1);
% AIC = AIC/AIC(1);
% AIC_c = AIC_c/AIC_c(1);
% cumulative = cumulative/cumulative(1);

% test = (2.*(1:p).*((1:p)+1))./(-(1:p)-1+N);

figure(), hold on;
p_axis = 1:p;
plot(p_axis, MDL)
plot(p_axis, AIC)
plot(p_axis, AIC_c)
plot(p_axis, cumulative)
legend('Minimum Description Length', 'Akaike Information Criterion', 'corrected AIC', 'Cumulative Error Squared')
title("Different criteria for loss calculation with respect to the model's order")
xlabel("model order p")
ylabel("Loss function")

figure(), hold on;
p_axis = 1:p;
plot(p_axis, MDL)
plot(p_axis, AIC)
plot(p_axis, AIC_c)
plot(p_axis, cumulative)
legend('Minimum Description Length', 'Akaike Information Criterion', 'corrected AIC', 'Cumulative Error Squared')
ylim([0.9 1])
title("Different criteria for loss calculation with respect to the model's order")
xlabel("model order p")
ylabel("Loss function")

%% Part 5
addpath(genpath([fileparts(pwd), '\Subpart1']));

for i = 1:p
    [h(i,:),w(i,:)]=freqz(1,[1 estimated_a(i,:)], N);
    per = pgm(estimated_signal(:,i));
    figure(), hold on
    plot(w(i,:)/(pi),per)
end

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
N = length(sun);

estimated_a = zeros(10);
aryule_coeff = zeros(10,11);

for p = 1:10 % order
    H = zeros(N-(p+1), p);
    for i = p+1:N
        for ii = 1:p
            H(i-p, ii) = sun(i-ii,1);
        end
    end
    
    estimated_a(p,1:p) = inv(H'*H)*H'*sun(p+1:end);
    aryule_coeff(p, 1:p+1) = aryule(sun, p); 
end

%% Part 4
N = length(sun);
sun_matrix = ones(N, p) .* sun;

estimated_signal = zeros(N, p);
for i = 1:p
    estimated_signal(:, i) = filter([0 estimated_a(i, :)], 1, sun);
end

e = zeros(1, p);

for i = 1:N
    e(1,:) = (estimated_signal(i,:) - sun_matrix(i,:)).^2 + e(1,:);
end

MDL = log(e(1,:)) + (1:p).*log(N)./N;
AIC = log(e(1,:)) + 2.*(1:p)./N;
AIC_c = AIC + (2.*(1:p).*((1:p)+1))./(-(1:p)-1+N);
cumulative = log(e(1,:));

% normalizing
MDL = MDL/MDL(1);
AIC = AIC/AIC(1);
AIC_c = AIC_c/AIC_c(1);
cumulative = cumulative/cumulative(1);

figure(), hold on;
p_axis = 1:p;
plot(p_axis, MDL)
plot(p_axis, AIC)
plot(p_axis, AIC_c)
plot(p_axis, cumulative)
legend('Minimum Description Length', 'Akaike Information Criterion', 'corrected AIC', 'Cumulative Error Squared')
title("Approximation error with respect to the model order")
xlabel("model order p")
ylabel("Loss function")

%% Part 5
addpath(genpath([fileparts(pwd), '\Subpart1']));
figure(), hold on
per = pgm(sun);
plot(linspace(0, 0.5, length(per)/2), per(1: end/2))
orders = [1 2 5 10];

for i = orders
    per = pgm(estimated_signal(:,i));
    plot(linspace(0, 0.5, length(per)/2), per(1:end/2))
end
title('PSD estimate of sunspot data')
xlabel('normalized frequency')
ylabel('Amplitude (Au)')
legend('Estimated periodogram for sunspot', 'Estimated periodogram for AR(1)', 'Estimated periodogram for AR(2)', 'Estimated periodogram for AR(5)', 'Estimated periodogram for AR(10)')
xlim([0 0.25])

figure(), hold on
for i = [2 9]
    per = pgm(estimated_signal(:,i));
    plot(linspace(0, 0.5, length(per)/2), per(1:end/2))
end
title('PSD estimate of sunspot data')
xlabel('normalized frequency')
ylabel('Amplitude (Au)')
legend('Estimated periodogram for AR(2)', 'Estimated periodogram for AR(9)')
xlim([0 0.25])

%% Part 6
N = 10:5:250;
mse = [];

for n = N
    H = zeros(n-3, 2);
    for i = 3:n
        for ii = 1:2
            H(i-2, ii) = sun(i-ii,1);
        end
    end
    
    a_est = inv(H'*H)*H'*sun(2+1:n);
    
    estimated_signal = filter([0; a_est], 1, sun);
    mse(find(N == n)) = mean((estimated_signal - sun).^2);
end

figure()
plot(N, mse)
title('MSE of the estimated series using different data lengths')
xlabel('data length (samples)')
ylabel('MSE')



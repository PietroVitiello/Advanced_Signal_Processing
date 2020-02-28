clear all
close all

load NASDAQ.mat
prices = NASDAQ.Close;
dates = NASDAQ.Date;

%% Part 1 pcf
p = 10;
pcf_axis = 1:p;

[a, ~, k] = aryule(prices, 10);
figure()
stem(pcf_axis, -k)

%% Part 1 errors
%% Part 4
N = length(prices);
p = 20;
norm_prices = zscore(prices);
close_p = ones(N, p) .* norm_prices;

filtered_prices = zeros(N, p);
for i = 1:p
    a = aryule(norm_prices, i);
    filtered_prices(:, i) = filter(-a, 1, norm_prices);
end

e = zeros(1, p);

for i = 1:N
    e(1,:) = (filtered_prices(i,:) - close_p(i,:)).^2 + e(1,:);
end

MDL = log(e(1,:)) + (1:p).*log(N)./N;
AIC = log(e(1,:)) + 2.*(1:p)./N;
AIC_c = AIC + (2.*(1:p).*((1:p)+1))./(-(1:p)-1+N);
% cumulative = log(e(1,:));

%normalizing
MDL = MDL/MDL(1);
AIC = AIC/AIC(1);
AIC_c = AIC_c/AIC_c(1);
% cumulative = cumulative/cumulative(1);

figure(), hold on;
p_axis = 1:p;
plot(p_axis, MDL)
plot(p_axis, AIC)
plot(p_axis, AIC_c)
% plot(p_axis, cumulative)
legend('Minimum Description Length', 'Akaike Information Criterion', 'corrected AIC') %, 'Cumulative Error Squared')
title("Different criteria for loss calculation with respect to the model's order")
xlabel("model order p")
ylabel("Loss function")
clear all
close all

load NASDAQ.mat
prices = zscore(NASDAQ.Close);
dates = NASDAQ.Date;
N = length(prices);

plot(prices)
a1 = aryule(prices, 1);

%% Part 1 pcf
p = 10;
pcf_axis = 1:p;

threshold = 1.96/sqrt(N);

[a, ~, k] = aryule(prices, 10);
figure(), hold on
stem(pcf_axis, -k)
plot([1 p], [1 1]'*[threshold -threshold], 'r')
title("Partial Correlation of NASDAQ prices data")
xlabel("\tau")
ylabel("PCF")

%% Part 1 errors
p = 10;
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
title("Loss calculation with respect to the model's order for NASDAQ data")
xlabel("model order p")
ylabel("Loss function")

%% Part 3
close all

sigma2 = 1:50:1001;
n = sigma2;

acf = xcorr(prices, 'unbiased');

for i = 1:length(sigma2)
    formula_sigma(i,:) = 1./n * (2*sigma2(1, i)^2);
end

figure()
h = heatmap(n, sigma2, formula_sigma);
h.Colormap = parula;
h.ColorScaling = 'log';
title("Heatmap for the variance of the noise variance")
xlabel("N")
ylabel("True variance of the driving noise")

for i = 1:length(sigma2)
    formula_a(i,:) = 1./n * (1/max(acf)) * sigma2(1, i);
end

figure()
h2 = heatmap(n, sigma2, formula_a);
h2.Colormap = parula;
h2.ColorScaling = 'log';
title('Heatmap for the variance of a1')
xlabel("N")
ylabel("True variance of the driving noise")



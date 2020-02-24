%% set up
clear all;
a1 = rand(100, 1).*5 - 2.5;
a2 = rand(100, 1).*3 - 1.5;
w = randn(1000, 1);

%% Part 1
clear pairs
tol = 10000;
pair_idx = 1;

for i = 1:length(a1)
    for j = 1:length(a2) 
        a = [1 -a1(i) -a2(j)];
        x = filter(1, a, w);
        
        if (x(end) < tol) && (x(end) > -tol)
            pairs(pair_idx,:) = [a1(i), a2(j)];
            pair_idx = pair_idx+1;
        end
    end 
end

scatter(pairs(:,1), pairs(:,2), 'r*')
title("Stable coefficients pairs for 100 coefficient values")
xlabel("a_{1}")
ylabel("a_{2}")
xlim([-2.5 2.5])
ylim([-1.5 1.5])

%% Part 2
close all
load sunspot.dat

lengths = [5 20 250];
index = 1;

for i = lengths
   sunspot_zero = sunspot(:,2) - mean(sunspot(1:i,2));
   acf = xcorr(sunspot(1:i, 2), 'unbiased');
   acf_zero = xcorr(sunspot_zero(1:i), 'coeff');
   figure()
   stem (-(i-1):(i-1), acf)
   title(sprintf("ACF of sunspot data for %.0f samples", i))
   ylabel("ACF")
   xlabel("Lag")
   figure()
   stem (-(i-1):(i-1), acf_zero)
   title(sprintf("ACF of zero-centered sunspot data for %.0f samples", i))
   ylabel("ACF")
   xlabel("Lag")
   series_means(index) = mean(sunspot(1:i, 2));
   index = index + 1;
end

%% DC offset
sunspot_zero = sunspot(1:20,2) - mean(sunspot(1:20,2));
const = ones(20, 1)*mean(sunspot(1:20,2));
const_acf = xcorr(const, 'unbiased');
acf_zero = xcorr(sunspot_zero, 'unbiased');
acf_complete = xcorr(sunspot_zero+const, 'unbiased');
figure()
stem(-19:19, const_acf)
figure()
stem(-19:19, const_acf+acf_zero)
figure()
stem(-19:19, acf_complete)

%% Part 3
close all;
clear all;

load sunspot.dat
sunspot_zero = sunspot(:,2) - mean(sunspot(:,2));

boundary = 1.96/sqrt(length(sunspot(:, 2)));

p = 10;
[a, ~, k] = aryule(sunspot(:,2), p);
[a_zero, ~, k_zero] = aryule(sunspot_zero, p);
pac = -k;
pac_zero = -k_zero;
figure(), hold on;
stem(1:p, pac)
stem(1:p, pac_zero)
plot([1 p], [1 1]' * [boundary -boundary], 'k--');
legend('With DC offset', 'Without DC offset')









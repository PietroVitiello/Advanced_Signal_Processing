%% Creation of signal
X = rand(1000, 1);  %Uniformly distributed between 0 and 1
x_axis = 1:1000;

figure()
plot(x_axis, X)
xlabel("discrete time (n)");

%% Part 1
sample_mean = mean(X);

%% Part 2
sample_standard_deviation = std(X);

%% Part 3
realizations = 1:10;
x_realizations = zeros(1000, 10);
sample_mean_realizations = zeros(1, 10);
std_realizations = zeros(1, 10);
bias_realizations = zeros(1, 10);

for realization = realizations
    x_realizations(:, realization) = rand(1, 1000);
    sample_mean_realizations(1, realization) = mean(x_realizations(:, realization));
    std_realizations(1, realization) = std(x_realizations(:, realization));
    bias_realizations(1, realization) = 0.5 - sample_mean_realizations(1, realization);
end

figure()
subplot(2,1,1)
scatter(realizations, sample_mean_realizations, 18,'ro', 'filled')
title("realizations' sample means")
ylabel("sample mean")
xlabel("realization")
ylim([0.45 0.55])
hline = refline([0, 0.5]);

subplot(2,1,2)
scatter(realizations, std_realizations, 18, 'ro', 'filled')
title("realizations' sample standard deviation")
ylabel("sample standard deviation")
xlabel("realization")
ylim([0.25 0.35])
hline = refline([0, 1/sqrt(12)]);

%% Part 4







